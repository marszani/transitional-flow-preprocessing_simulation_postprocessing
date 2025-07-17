#!/bin/bash

CONTAINER="mach_aero_backup"
PYTHON_BIN="/home/mdolabuser/.pyenv/versions/3.9.12/bin/python"
SCRIPT_PATTERN=".*\.py"

# Determine signal based on input argument
case "$1" in
  save)
    SIGNAL="USR1"           # save simulation (your current signal for save)
    ;;
  stop)
    SIGNAL="TERM"           # graceful stop
    ;;
  kill)
    SIGNAL="KILL"           # force kill
    ;;
  save_and_stop)
    SIGNAL="USR2"           # save and stop simulation
    ;;
  *)
    echo "Usage: $0 {save|stop|kill|save_and_stop}"
    exit 1
    ;;
esac

echo "Looking for processes with python interpreter: $PYTHON_BIN"
echo "Matching script pattern: $SCRIPT_PATTERN"
echo "Sending signal: SIG$SIGNAL"

# Get matching PIDs: command must be exact python interpreter, and script must match pattern
PIDS=$(docker exec "$CONTAINER" ps aux | awk -v py="$PYTHON_BIN" -v sp="$SCRIPT_PATTERN" '$11 == py && $0 ~ sp { print $2 }')

if [ -z "$PIDS" ]; then
  echo "No matching processes found."
  exit 0
fi

# Send the signal to each PID
for pid in $PIDS; do
  echo "Sending SIG$SIGNAL to PID $pid"
  docker exec "$CONTAINER" kill -"$SIGNAL" "$pid"
done