import subprocess

subprocess.check_call(["pip", "install", "-r", "requirements.txt"])

with open("packages_loaded.txt", "w") as f:
    f.write("package loading is done!")
