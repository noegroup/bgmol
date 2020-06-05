
import subprocess as sp
import os


def download_file(host, path, destination, protocol="rsync"):
    command = f"rsync --info=progress2 -a {host}:{path} {destination}"
    print(command)
    os.system(command)


