"""
Graph-based genetic algorithm (after Jensen 2019, "GB-GA").

A genetic algorithm imitates selection: keep a population of molecules, let the
fitter ones "breed" (crossover) and "mutate" to propose children, score the
children, keep the best. Repeat. There is no gradient and no learning; the only
thing steering the search is the fitness number, which here comes from the
surrogate. That is why it is a clean test of whether the surrogate can be
exploited.

One generation:
  1. Pick two parents, favouring high fitness (roulette-wheel selection).
  2. Child = crossover of the parents; with probability GA_MUTATION_PROB, or if
     crossover failed, the child is mutated as well.
  3. The child is prepared under the chosen policy (candidates.py), deduplicated
     against the current population and this batch, and scored.
  4. Merge parents and children, keep the best GA_POP_SIZE by scalar fitness.
     Elitist: a molecule with a top fitness is never lost. That guarantees retention by fitness, NOT
     that every component of the objective improves monotonically (the composite can rise while
     one term falls). This is greedy on purpose, so exploitable surrogate weaknesses get exploited.

"Unseen" means only: not in the current population and not already proposed in this batch. A proposal can
have been generated or scored in an earlier generation, in another run, or be a training molecule, so it
is not a claim of global novelty. No global novelty filter is applied, since that would change the search.

Similarity floor (H5). With min_similarity set, a molecule whose max-Tanimoto to the training set is below
the floor (or not finite) is INFEASIBLE: it is removed before survivors are chosen, whatever its fitness. This
is distinct from a feasible molecule whose fitness happens to be zero, which can still survive.
"""

import math

import numpy as np
import pandas as pd
from rdkit import Chem

from cdk2moo import config
from cdk2moo.candidates import POLICIES, prepare_candidate
from cdk2moo.mutations import crossover, mutate
from cdk2moo.objectives import score_molecules, fitness


def _count(rejects, key, n=1):
    rejects[key] = rejects.get(key, 0) + n


def _check_floor(min_similarity):
    """A similarity floor is None (off) or a finite number in [0, 1] (the range of a Tanimoto)."""
    if min_similarity is None:
        return
    if isinstance(min_similarity, bool) or not isinstance(min_similarity, (int, float)) \
            or not math.isfinite(min_similarity) or not 0.0 <= min_similarity <= 1.0:
        raise ValueError(f"min_similarity must be None or a finite number in [0, 1], got {min_similarity!r}")


def _feasible(scores, min_similarity):
    """Boolean mask: max-Tanimoto is finite and at least the floor. All True when the floor is off."""
    if min_similarity is None:
        return pd.Series(True, index=scores.index)
    similarity = scores["max_tanimoto"].astype(float)
    return np.isfinite(similarity) & (similarity >= min_similarity)


def _propose_children(population, known, rng, size_range, policy="legacy", cache=None, rejects=None):
    """
    Generate up to GA_N_CHILDREN new molecules that are valid under the policy and not already known.

    Returns a list of dicts: smiles (the representation that will be scored), proposal_smiles (the edited
    graph before preparation), how it was made ("crossover", "crossover+mutation" or "mutation") and its
    parents' SMILES (parent_b is empty for a pure mutation).
    """
    closed_shell = policy == "corrected"
    mols = [Chem.MolFromSmiles(s) for s in population["smiles"]]
    # Fitness-proportional selection; a tiny floor keeps every molecule eligible.
    weights = population["fitness"].to_numpy() + 1e-6
    weights = weights / weights.sum()

    children = []
    seen = set(known)
    attempts = 0
    while len(children) < config.GA_N_CHILDREN and attempts < 50 * config.GA_N_CHILDREN:
        attempts += 1
        a, b = rng.choice(len(mols), size=2, replace=False, p=weights)
        child = crossover(mols[a], mols[b], rng, size_range, closed_shell=closed_shell, rejects=rejects)
        origin = "crossover"
        if child is None or rng.random() < config.GA_MUTATION_PROB:
            origin = "mutation" if child is None else "crossover+mutation"
            child = mutate(child if child is not None else mols[a], rng, size_range,
                           closed_shell=closed_shell, rejects=rejects)
        if child is None:
            continue
        _count(rejects, "proposals_attempted")
        smiles, reason, status = prepare_candidate(child, policy, size_range, cache=cache, rejects=rejects)
        if smiles is None or smiles in seen:
            if smiles is not None and rejects is not None:
                _count(rejects, "duplicate")
            continue
        seen.add(smiles)
        _count(rejects, "proposals_accepted")
        children.append({
            "smiles": smiles,
            "tautomer_status": status or "not_applicable",
            "proposal_smiles": Chem.MolToSmiles(child),
            "origin": origin,
            "parent_a": population["smiles"].iloc[a],
            "parent_b": population["smiles"].iloc[b] if origin != "mutation" else "",
        })
    return children


def _prepare_start(start_smiles, policy, size_range, cache, rejects):
    """
    Start molecules as the representation that will be scored.

    legacy: canonical SMILES, as in every historical run. corrected: standardized and validated like any
    other candidate; the population must not shrink, so any failure or collision raises instead of
    silently starting smaller than requested.
    """
    mols = [Chem.MolFromSmiles(s) for s in start_smiles]
    if policy == "legacy":
        return [Chem.MolToSmiles(m) for m in mols]
    prepared, problems = [], []
    for original, mol in zip(start_smiles, mols):
        smiles, reason, _ = prepare_candidate(mol, policy, size_range, cache=cache, rejects=rejects)
        if smiles is None:
            problems.append(f"{original}: {reason}")
        prepared.append(smiles)
    if problems:
        raise ValueError(f"{len(problems)} of {len(prepared)} starting molecules fail the {policy!r} policy, e.g. {problems[:3]}")
    if len(set(prepared)) < len(prepared):
        raise ValueError(f"standardization merged {len(prepared) - len(set(prepared))} starting molecules into duplicates; "
                         "the initial population would be smaller than requested")
    return prepared


def run_ga(start_smiles, arm, real_forest, scrambled_forest, train_fps, seed,
           size_range, n_generations=None, min_similarity=None, policy="legacy", rejects=None):
    """
    Run one GA and return every generation's population as one DataFrame
    (columns: the scores from objectives.score_molecules, plus fitness, origin, proposal_smiles, tautomer_status,
    parent_a, parent_b, birth_generation, generation, arm, seed, min_similarity, policy). Every molecule keeps its proposal, scored
    string, parents, origin, birth generation, seed, policy and the tautomer-enumeration status of its standardization; rejection reasons and
    the numbers of proposals attempted and accepted are in history.attrs["rejects"].

    The forests are passed in already trained and are never modified here.
    `seed` controls only the search (parent choice, edits), so two arms given
    the same seed and start molecules differ only in their fitness function.

    policy: "legacy" (default, every historical run) or "corrected" (standardize and validate every candidate
    like the training molecules, closed-shell only); see candidates.py. min_similarity: see the module text; the
    starting molecules must already satisfy it or a ValueError is raised (the floor is never silently relaxed).
    `rejects`, a dict, accumulates rejection counts by reason; the same dict is available afterwards in
    `history.attrs["rejects"]`.
    """
    if policy not in POLICIES:
        raise ValueError(f"unknown preparation policy {policy!r}; choose from {POLICIES}")
    _check_floor(min_similarity)
    if n_generations is None:
        n_generations = config.GA_N_GENERATIONS
    rejects = {} if rejects is None else rejects
    cache = {}
    rng = np.random.default_rng(seed)

    start = _prepare_start(list(start_smiles), policy, size_range, cache, rejects)
    population = score_molecules(start, real_forest, scrambled_forest, train_fps)
    population["fitness"] = fitness(population, arm)
    feasible = _feasible(population, min_similarity)
    if not feasible.all():
        raise ValueError(f"{int((~feasible).sum())} of {len(population)} starting molecules violate the similarity floor "
                         f"{min_similarity}; refusing to start (the floor is not relaxed and infeasible molecules are not used)")
    population["origin"] = "start"
    population["tautomer_status"] = [prepare_candidate(Chem.MolFromSmiles(x), policy, size_range, cache=cache)[2] or "not_applicable"
                                     for x in population["smiles"]]
    population["proposal_smiles"] = population["smiles"]
    population["parent_a"] = ""
    population["parent_b"] = ""
    population["birth_generation"] = 0
    history = [population.assign(generation=0)]

    for generation in range(1, n_generations + 1):
        children = _propose_children(population, population["smiles"], rng, size_range, policy, cache, rejects)
        if children:
            info = pd.DataFrame(children)
            scored = score_molecules(list(info["smiles"]), real_forest, scrambled_forest, train_fps)
            scored["fitness"] = fitness(scored, arm)
            scored["origin"] = info["origin"]
            scored["tautomer_status"] = info["tautomer_status"]
            scored["proposal_smiles"] = info["proposal_smiles"]
            scored["parent_a"] = info["parent_a"]
            scored["parent_b"] = info["parent_b"]
            scored["birth_generation"] = generation
            feasible_children = _feasible(scored, min_similarity)
            if not feasible_children.all():
                _count(rejects, "infeasible_similarity", int((~feasible_children).sum()))
            pool = pd.concat([population, scored[feasible_children]], ignore_index=True)
        else:
            pool = population
        if len(pool) < config.GA_POP_SIZE:
            raise RuntimeError(f"only {len(pool)} feasible molecules for a population of {config.GA_POP_SIZE}")
        population = pool.sort_values("fitness", ascending=False).head(config.GA_POP_SIZE)
        population = population.reset_index(drop=True)
        history.append(population.assign(generation=generation))

    out = pd.concat(history, ignore_index=True).assign(arm=arm, seed=seed, min_similarity=min_similarity, policy=policy)
    out.attrs["rejects"] = dict(rejects)
    return out
