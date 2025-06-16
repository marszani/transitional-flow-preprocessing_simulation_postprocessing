# Simulation Control Script for Docker Container

**File:** `run_all_sims.sh`

## Purpose

Control running Python simulation processes inside a Docker container by sending signals for saving, stopping, or killing the simulations gracefully or forcibly.

---

## Usage

```bash
./run_all_sims.sh {save|stop|kill|save_and_stop}

	•	save: Sends a user-defined signal (USR1) to request simulation save without stopping.
	•	stop: Sends a termination signal (TERM) for graceful shutdown.
	•	kill: Sends a kill signal (KILL) to forcibly terminate simulations.
	•	save_and_stop: Sends a user-defined signal (USR2) to request save and then stop.

⸻

Details
	•	Container Name:
The script targets the Docker container named "mach_aero_backup".
	•	Python Interpreter Path:
It looks for processes running the exact Python interpreter at:
/home/mdolabuser/.pyenv/versions/3.9.12/bin/python.
	•	Script Pattern:
Matches running Python scripts whose command line matches the regex pattern .*\.py.
	•	Process Identification:
Uses ps aux inside the container combined with awk to identify PIDs of Python processes matching the pattern.
	•	Signal Sending:
Sends the appropriate UNIX signal to each matched PID inside the Docker container.

⸻

Notes
	•	Requires Docker installed and the specified container running.
	•	The Python interpreter path inside the container must match the one in the script.
	•	The user running the script must have permission to execute Docker commands.
	•	The script prints status messages indicating the processes found and signals sent.

⸻

Example

To gracefully stop all running simulations and allow cleanup:

./run_all_sims.sh stop

To force kill all simulations immediately:

./run_all_sims.sh kill

