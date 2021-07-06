
import subprocess
import psutil

def run_output(command):
    print(command)

    command_split = list(map(lambda x:x.replace("#", " "),command.split(" ")))
    proc = subprocess.Popen(command_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    _stdout, _stderr = proc.communicate()
    stdout = _stdout.decode('utf-8').replace("\\n","\n")
    stderr = _stderr.decode('utf-8').replace("\\n","\n")
    returncode = proc.returncode
    

    print("Current number of open files: ",len(psutil.Process().open_files()))
    return stdout, stderr, returncode
