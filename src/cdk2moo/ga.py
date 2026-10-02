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
  3. Score the new children.
  4. Merge parents and children, drop duplicates, keep the best GA_POP_SIZE.
     Elitist: a good molecule never disappears. This is greedy on purpose, so
     any exploitable weakness in the surrogate gets exploited.
"""

import numpy as np
import pandas as pd
from rdkit import Chem

from cdk2moo import config
from cdk2moo.mutations import crossover, mutate
from cdk2moo.objectives import score_molecules, fitness


def _propose_children(population, known, rng, size_range):
    """
    Generate GA_N_CHILDREN new, valid, previously unseen molecules.

    Returns a list of dicts: the child's SMILES, how it was made ("crossover",
    "crossover+mutation" or "mutation") and its parents' SMILES (parent_b is
    empty for a pure mutation).
    """
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
        child = crossover(mols[a], mols[b], rng, size_range)
        origin = "crossover"
        if child is None or rng.random() < config.GA_MUTATION_PROB:
            origin = "mutation" if child is None else "crossover+mutation"
            child = mutate(child if child is not None else mols[a], rng, size_range)
        if child is None:
            continue
        smiles = Chem.MolToSmiles(child)
        if smiles in seen:
            continue
        seen.add(smiles)
        children.append({
            "smiles": smiles,
            "origin": origin,
            "parent_a": population["smiles"].iloc[a],
            "parent_b": population["smiles"].iloc[b] if origin != "mutation" else "",
        })
    return children


def run_ga(start_smiles, arm, real_forest, scrambled_forest, train_fps, seed,
           size_range, n_generations=None, min_similarity=None):
    """
    Run one GA and return every generation's population as one DataFrame
    (columns: the scores from objectives.score_molecules, plus fitness,
    origin, parent_a, parent_b, birth_generation, generation, arm, seed).

    The forests are passed in already trained and are never modified here.
    `seed` controls only the search (parent choice, edits), so two arms given
    the same seed and start molecules differ only in their fitness function.

    `min_similarity` (default None = off) is a hard floor for the H5 mitigation arm: a
    molecule whose max-Tanimoto to the training set is below it gets fitness 0 and so cannot
    survive selection.
    """
    if n_generations is None:
        n_generations = config.GA_N_GENERATIONS
    rng = np.random.default_rng(seed)
    start = [Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in start_smiles]
    population = score_molecules(start, real_forest, scrambled_forest, train_fps)
    population["fitness"] = fitness(population, arm)
    if min_similarity is not None:
        population.loc[population["max_tanimoto"] < min_similarity, "fitness"] = 0.0
    population["origin"] = "start"
    population["parent_a"] = ""
    population["parent_b"] = ""
    population["birth_generation"] = 0
    history = [population.assign(generation=0)]

    for generation in range(1, n_generations + 1):
        children = _propose_children(population, population["smiles"], rng, size_range)
        if children:
            info = pd.DataFrame(children)
            scored = score_molecules(list(info["smiles"]), real_forest,
                                     scrambled_forest, train_fps)
            scored["fitness"] = fitness(scored, arm)
            if min_similarity is not None:
                scored.loc[scored["max_tanimoto"] < min_similarity, "fitness"] = 0.0
            scored["origin"] = info["origin"]
            scored["parent_a"] = info["parent_a"]
            scored["parent_b"] = info["parent_b"]
            scored["birth_generation"] = generation
            pool = pd.concat([population, scored], ignore_index=True)
        else:
            pool = population
        population = pool.sort_values("fitness", ascending=False).head(config.GA_POP_SIZE)
        population = population.reset_index(drop=True)
        history.append(population.assign(generation=generation))

    return pd.concat(history, ignore_index=True).assign(arm=arm, seed=seed, min_similarity=min_similarity)
