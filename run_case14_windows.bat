@echo off
setlocal
cd /d %~dp0
if not exist .venv (
  python -m venv .venv
)
call .venv\Scripts\activate
python -m pip install -r requirements.txt
python run_experiments.py --cases case14 --repeats 3 --tol 1e-8
python export_latex_tables.py
pause
