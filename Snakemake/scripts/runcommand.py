
import subprocess


def run_output(command):
    print(command)
    command_split = list(map(lambda x:x.replace("#", " "),command.split(" ")))
    proc = subprocess.Popen(command_split, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode('utf-8').replace("\\n","\n")
    stderr = stderr.decode('utf-8').replace("\\n","\n")
    returncode = proc.returncode
    
    return stdout, stderr, returncode
