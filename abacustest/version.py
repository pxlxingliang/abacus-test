import subprocess
import os

BASE_VERSION = "v0.4.50"

def _no_need_subversion():
    try:
        modified_files = subprocess.check_output(
            ["git", "diff-tree", "--no-commit-id", "--name-only", "-r", "HEAD"],
            cwd=os.path.dirname(os.path.abspath(__file__)),
            stderr=subprocess.DEVNULL
        ).decode().splitlines()
        
        current_file = os.path.relpath(__file__, os.getcwd())
        
        # If this file is modified, indicating the main version is updated
        if current_file in modified_files:
            return True
        
        # If codes in the abacustest directory are modified, indicating the main version is not updated
        for i in modified_files:
            if i.startswith("abacustest/"):
                return False
        
        return True
    except:
        return False

def get_version():
    try:
        if _no_need_subversion():
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