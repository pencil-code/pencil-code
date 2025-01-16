#!/usr/bin/env python
# coding: utf-8
import pencil as pc
from pencil.math import natural_sort
from pencil.io import mkdir
import os
from os.path import join
import glob
import numpy as np
import time
import inspect
import odop
import cv2
import glob


def get_reduce_task(field_name):
    @odop.task(name=f"reduce_{field_name}",trigger="file_updated",file_path="/users/toukopur/pencil-code/samples/gputest/data/proc0/var.dat")
    def reduce(filenames=None):
        f = pc.read.var()
        #file_heavy = False
        field = getattr(f,field_name)
        for i in range(100):
            xy_aver = np.sum(field,2)
            np.save(f"/users/toukopur/odop-dst/{field_name}_xyaver.dat",xy_aver)

        print(f"{field_name}_xyaver: DONE\n");
    return reduce

for field_name in [
                    "ux",
                    "uy",
                    "uz",
                    "ax",
                    "ay",
                    "az",
                    "lnrho",
                    "ss"
                    ]:
        get_reduce_task(field_name)



import os
import swiftclient
from swiftclient.multithreading import OutputManager
from swiftclient.service import SwiftError, SwiftService, SwiftUploadObject


def get_conn():
    _authurl = os.environ['OS_AUTH_URL']
    _auth_version = os.environ['OS_IDENTITY_API_VERSION']
    _user = os.environ['OS_USERNAME']
    _key = os.environ['OS_PASSWORD']
    _os_options = {
        'user_domain_name': os.environ['OS_USER_DOMAIN_NAME'],
        'project_domain_name': os.environ['OS_USER_DOMAIN_NAME'],
        'project_name': os.environ['OS_PROJECT_NAME']
    }
    
    return swiftclient.Connection(
        authurl=_authurl,
        user=_user,
        key=_key,
        os_options=_os_options,
        auth_version=_auth_version
    )
container = "NEOSC"

def download():
    conn = get_conn()
    resp_headers, containers = conn.get_account()
    my_obj = conn.get_container(container, prefix="allas/")[1]
    for obj in my_obj:
        object_name = obj["name"]
        # Define the local path by removing the root prefix from the object name
        local_path = os.path.join(object_name)

        # Create the necessary directories for nested folders
        os.makedirs(os.path.dirname(local_path), exist_ok=True)

        # Download and save each object to its corresponding local path
        with open(local_path, 'wb') as out_file:
            headers, object_contents = conn.get_object(container, object_name)
            out_file.write(object_contents)
            print(f'Downloaded {object_name} to {local_path}')

def upload():
    conn = get_conn()
    resp_headers, containers = conn.get_account()
    src_name = "output.mp4"
    dst_name = "output.mp4"
    container = "NEOSC"
    with open(src_name,'rb') as file_data:
        conn.put_object(container,dst_name, contents = file_data, content_type ='video/mp4')
    print(f"Uploaded {dst_name} to container {container} as {src_name}")

def create_video(image_folder, output_file, framerate=30):
    images = sorted(glob.glob(f"{image_folder}/*.png"))
    frame = cv2.imread(images[0])
    height, width, layers = frame.shape

    video = cv2.VideoWriter(output_file, cv2.VideoWriter_fourcc(*'mp4v'), framerate, (width, height))
    for image in images:
        video.write(cv2.imread(image))

    video.release()
    print(f"Video saved as {output_file}")

@odop.task(name="download_and_render",max_runs=1)
def main():
    download()
    datadir="/users/toukopur/pencil-code/samples/gputest/allas/data"
    import pencil as pc
    slices = pc.read.slices(datadir=datadir)
    fields = ["uu1","uu2"]
    pc.visu.rvid_box.plot_box(slices,datadir=datadir,fields=fields,colorscale="Jet")
    create_video("images", "output.mp4", framerate=1)
    upload()

if __name__ == "__main__":
    main()
    #upload()
