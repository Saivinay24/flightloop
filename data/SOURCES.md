# Data sources (all public or synthetic)

| path | source | URL | page / table | method | date | license |
|---|---|---|---|---|---|---|
| tm4062/tm4062.pdf | NASA TM-4062, McGhee, Walker, Millard 1988, public-experimental | https://ntrs.nasa.gov/api/citations/19890001471/downloads/19890001471.pdf | whole report | download | 2026-11-18 | US Government work, public domain |
| tm4062/tm4062.txt | text layer of tm4062.pdf | derived | all | pdftotext -layout | 2026-11-18 | as above |
| tm4062/pass_b/polar_re200k.csv | TM-4062 Table B1, runs 9,10,13, R = 200,000 | as above | doc p. 92 (PDF p. 96) | manual transcription from 150 dpi page image (pass B) | 2026-11-18 | as above |
| tm4062/pass_b/polar_re100k.csv | TM-4062 Table B1, runs 15,16, R = 100,000 | as above | doc p. 93 (PDF p. 97) | manual transcription (pass B) | 2026-11-18 | as above |
| tm4062/pass_b/oilflow.csv | TM-4062 Table III oil flow LS/TR | as above | doc p. 11 (PDF p. 15) | manual transcription (pass B) | 2026-11-18 | as above |
| tm4062/pass_b/cp_upper_re200k.csv | TM-4062 Appendix D upper-surface Cp, runs 9,10,13, alpha 0.01/3.99/8.00 | as above | doc p. 126, 128, 130 (PDF p. 130, 132, 134) | manual transcription (pass B) | 2026-11-18 | as above |
| airfoils/* | UIUC Airfoil Coordinates Database (Selig), public-experimental geometry, copied from StratoGlide `data/validation/coords/` (its MANIFEST.csv has the per-file URLs) | https://m-selig.ae.illinois.edu/ads/coord_database.html | per airfoil .dat | copy | 2026-11-18 | UIUC database terms (free for research use) |
| tm4062/pass_a/* | TM-4062 pass A transcription (text layer + second manual read) | as TM-4062 | B1 PDF p. 96-97, App. D PDF p. 130-134, Table III PDF p. 15 | tools/ocr_tm4062.py + manual crops | 2026-11-18 | public domain |
| tm4062/canonical/* | resolved canonical TM-4062 tables used by all scripts | derived | as above | scripts/diff_transcriptions.py + RESOLUTIONS.csv | 2026-11-18 | public domain |
| tm4062/RESOLUTIONS.csv | page-image resolutions of every pass A / pass B disagreement | derived | as above | manual, 150/300 dpi page images | 2026-11-18 | public domain |
| tm4062/diff_report.json | output of scripts/diff_transcriptions.py | derived | | script | 2026-11-18 | |
| uiuc_lsat/* | UIUC Low-Speed Airfoil Tests, Summary of Low-Speed Airfoil Data vol. 1 (Selig, Guglielmo, Broeren, Giguere 1995) and vol. 3 (Selig, Lyon, Giguere, Ninham, Guglielmo 1996), public experimental; copied from StratoGlide `data/validation/uiuc_lsat/` (its MANIFEST.csv has SHA-256 per file) | https://m-selig.ae.illinois.edu/pd/pub/lsat/volume01.zip, volume03.zip | DRAG/LIFT/COORD text files | copy | 2026-11-18 | GPL (GPL.TXT shipped with the data) |
| alfa/* | CMU ALFA dataset (Keipour, Mousaei, Scherer), public experimental, CC0: processed.zip (CSV topics per sequence), dataflash.zip (ArduPlane .bin), README.txt; raw archives and extracted files are not committed (`.gitignore`), AIRCRAFT.md and avl/ are FlightLoop notes | https://doi.org/10.1184/R1/12707963 (figshare files 24095870, 24096047, 24098639) | engine-failure sequences | download | 2026-11-18 | CC0 |
| px4/* | PX4 Flight Review public logs 417157ed-a540-4e7d-9c85-a20fa5023e2b, 61bd5632-ee50-44af-a3c9-7ae2ce28865b, 4af54f87-a19b-4254-80b6-31e64617ff10 and the public index dbinfo.json; used for parser verification only (logs not committed) | https://review.px4.io/dbinfo, https://cdn.logs.px4.io/<id>.ulg | whole logs | download | 2026-11-18 | public Flight Review uploads |
