This folder's subfolders named PET and MR contain Python demo scripts for PET and MR reconstruction respectively.

For a brief description of each script's purpose and usage type
```
python <script_name>.py --help
```
(on VM use `python3` instead of `python`).

Note that the SIRF Python utilities should be in your `PYTHONPATH` (see the installation instructions).

To run the demo scripts, please create an environment variable `SIRF_PATH` whose value is the location of your SIRF source, e.g. `/home/sirfuser/devel/SIRF`. (NOTE: even on Windows, you must use `/` in paths, not `\\`.) If you do not do this, you can still run the demos but you will have to give full path to raw data files via command-line options `-p` or `--path`.

