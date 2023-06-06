import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging;
import time
import psutil
from datetime import datetime

logging.basicConfig(level=logging.DEBUG)

# Read the list of commands from a file
filename = "/tmp/executeTheseFastOutputs.txt"  # Update with your file name
with open(filename, "r") as file:
    commands = file.read().splitlines()

commandsToDo = commands
print(len(commandsToDo))

logging.debug(f"{str(datetime.now())}: Beginning ALL with RAM at {str(psutil.virtual_memory()[2])}")

# Function to execute a command and print the execution time
def execute_command(command):
    start_time = time.time()
    printcommand = command.split(" /tmp/")[0]
    logging.debug(f"{str(datetime.now())}: START'{printcommand}' with RAM at {str(psutil.virtual_memory()[2])}")
    process = subprocess.Popen(command, shell=True)
    process.communicate()  # Wait for the command to finish
    end_time = time.time()
    execution_time = end_time - start_time
    logging.debug(f"{str(datetime.now())}: DONE '{printcommand}' executed in {execution_time:.2f} seconds.")

# Execute commands in parallel using 12 cores
with ThreadPoolExecutor(max_workers=14) as executor:
    # Submit commands for execution
    futures = [executor.submit(execute_command, command) for command in commandsToDo]

    # Wait for all commands to complete
    for future in as_completed(futures):
        future.result()

print("Done.")
