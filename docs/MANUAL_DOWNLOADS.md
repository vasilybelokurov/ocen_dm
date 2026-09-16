# Manual data retrieval

Two of the four primary datasets are journal supplementary material and cannot
be fetched programmatically (publisher sessions / anti-bot pages). Retrieve them
by hand, put them where this page says, then run `ocen preprocess`. The loaders
verify the columns; nothing is filled in on your behalf.

`ocen fetch-data` reports these as `manual` rather than failing.

---

## 1. Kuzma et al. (2026) spectroscopy

- Paper: https://arxiv.org/abs/2605.23474
- DOI: https://doi.org/10.1093/mnras/stag1147
- Publisher: https://academic.oup.com/mnras/article/550/2/stag1147/8709286
- Underlying ESO data: programme `108.22MM.001`, PI Kuzma (not re-reduced in v1).

Steps:

1. Open the publisher page and download the online supplementary star list.
2. Save it as `data/raw/kuzma2026_spectroscopy/kuzma2026_members.<ext>`, keeping
   the original format (`.fits`, `.csv`, `.txt`, `.vot` are all readable).

   **If the file is plain CSV or ASCII it carries no units.** The loader will
   refuse to assume them -- it prints the columns and stops. Declare the input
   units in `configs/column_maps.yaml`, taking them from the paper's table
   caption, e.g.

   ```yaml
   kuzma2026_spectroscopy:
     vlos:       {column: HRV,  unit: km / s}
     vlos_error: {column: e_RV, unit: km / s}
   ```

   Record in `JOURNAL.md` which table caption you took each unit from.
3. Record the download date and the exact supplementary file name in
   `JOURNAL.md`.
4. Run:

   ```bash
   source ~/Work/venvs/.venv/bin/activate
   PYTHONPATH=src python -m ocen_dm.cli preprocess --dataset kuzma2026
   ```

If a required role (`ra`, `dec`, `vlos`, `vlos_error`, `source_id`) does not
resolve, the loader prints the table's real column names. Add the mapping to
`configs/column_maps.yaml` under `kuzma2026_spectroscopy`.

**Likelihood rule.** The sample is target-selected: use velocities and
metallicities conditional on the observed target positions. Spectroscopic
detection counts are not an unbiased tail surface density.

---

## 2. Fimbulthul (Ibata et al. 2019) -- no longer manual

- DOI: https://doi.org/10.1038/s41550-019-0751-x
- Publisher: https://www.nature.com/articles/s41550-019-0751-x

**This download is now automatic.** CDS/VizieR hosts the candidate list as
`J/other/NatAs/3.667/tables1` ("Candidate members of the Fimbulthul stellar
stream", 309 rows), so `ocen fetch-data` retrieves it directly:

```bash
PYTHONPATH=src python -m ocen_dm.cli fetch-data --dataset fimbulthul_ibata2019
PYTHONPATH=src python -m ocen_dm.cli preprocess --dataset fimbulthul
```

Two provenance caveats are recorded in `provenance/datasets.yaml`:

- it is a **CDS-curated representation** of Supplementary Table 1, not the
  publisher's own file: CDS renames the columns, adds `SimbadName` and `recno`;
- the VOTable header embeds the export timestamp, so a re-export has a
  different checksum. Compare row count and values, not bytes.

If you prefer the publisher's own file, download Supplementary Table 1 and save
it as `data/raw/fimbulthul_ibata2019/fimbulthul_supp_table1.<ext>` (only one
file may match that pattern), converting a PDF/xlsx to CSV **without editing
values** and noting the conversion command in `JOURNAL.md`.

**Likelihood rule.** STREAMFINDER selection is complex: use Fimbulthul as a
stream-track constraint only, not as an absolute surface density.

---

## 3. Zenodo datasets when the API is down

`zenodo.org` returned HTTP 504 for both the API and the file endpoints on
2026-09-16. `ocen fetch-data` retries with exponential backoff and reports
failures; it never writes a placeholder file. Re-run it when Zenodo is up:

```bash
PYTHONPATH=src python -m ocen_dm.cli fetch-data
```

If you must download by hand:

- Kuzma & Ishigaki (2025): https://zenodo.org/records/14791430 ->
  `data/raw/kuzma2025_pristine/wCen_table.fits`
  (published MD5 `3dcd58e69767901fbbdf69385a14c768`)
- oMEGACat VI: https://zenodo.org/records/14978551 ->
  `data/raw/omegacat_vi_kinematics/`

Then re-run `ocen fetch-data`, which verifies checksums of files already on
disk and writes `provenance/manifest.json` and `provenance/checksums.txt`.
