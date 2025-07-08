import subprocess
import os

BASE_VERSION = "v0.4.23"

def _is_version_file_updated():
    try:
        modified_files = subprocess.check_output(
            ["git", "diff-tree", "--no-commit-id", "--name-only", "-r", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            stderr=subprocess.DEVNULL
        ).decode().splitlines()
        
        current_file = os.path.relpath(__file__, os.getcwd())
        return current_file in modified_files
    except:
        return False

def get_version():
    try:
        if _is_version_file_updated():
            return BASE_VERSION
        
        commit_hash = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.STDOUT,
            text=True
        ).strip()
        

        timestamp = subprocess.check_output(
            ["git", "log", "-1", "--format=%cd", "--date=format:%Y%m%d%H%M%S"],
            stderr=subprocess.STDOUT,
            text=True
        ).strip()
        
        return f"{BASE_VERSION}+{timestamp}.{commit_hash}"
    except:
        return BASE_VERSION  

__version__ = get_version()