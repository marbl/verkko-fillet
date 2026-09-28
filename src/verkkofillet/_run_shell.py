import sys
import shlex
import subprocess
import os
import re
import shutil
from datetime import datetime

def run_shell(cmd: str, wkDir: str, functionName: str, longLog: bool = False, showOnly: bool = False) -> None:
    """
    Run a shell command and log its output to a file.

    Parameters:
        cmd (str): The shell command to execute.
        wkDir (str): Working directory for the command.
        functionName (str): Name of the function or task for logging.
        longLog (bool): Whether to include the log file path in success message.
        showOnly (bool): Only print cmd without running.

    Raises:
        subprocess.CalledProcessError: If the command execution fails.
    """
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(wkDir, "log", f"{functionName}.{now}.log")
    
    # Ensure the log directory exists
    os.makedirs(os.path.dirname(log_file), exist_ok=True)

    if showOnly : 
        print(cmd)
        return
        
    try:
        with open(log_file, "w") as log:
            log.write(f"Running command: {cmd}\n")
            subprocess.run(
                cmd,
                stdout=log,  # Redirect stdout to log file
                stderr=subprocess.STDOUT,  # Merge stderr with stdout
                shell=True,
                check=True,
                cwd=wkDir  # Set working directory
            )
        # Print success message based on `longLog` flag
        message = (
            f"[{functionName}] Command executed successfully. Logs saved to {log_file}."
            if longLog
            else f"[{functionName}] Command executed successfully!"
        )
        print(message)

    except subprocess.CalledProcessError as e:
        # Enhanced error reporting for subprocess failures
        print(f"[{functionName}] Command failed: {cmd}")
        print(f"[{functionName}] Error code: {e.returncode}")
        print(f"[{functionName}] Check logs for details: {log_file}.")
        raise
    
    except Exception as e:
        # Generic error handling for other unexpected exceptions
        print(f"[{functionName}] An unexpected error occurred: {e}")
        raise