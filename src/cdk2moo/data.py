"""Stage 1a: fetching bioactivity data from ChEMBL.

Will live here:

- resolution of the CDK2 target ID through the ChEMBL web client, checked
  against ``config.EXPECTED_CHEMBL_TARGET_ID`` rather than assumed
- download of activity records (IC50 / Ki) for that target
- caching of every raw response under ``config.RAW_DIR`` so that re-runs do
  not hit the API, which is slow

Raw responses are written verbatim. All cleaning happens in ``standardize``.
"""
