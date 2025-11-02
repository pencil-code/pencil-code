.. _howtodocker:

*********************************
Super short Docker reference
*********************************

In this age of instant internet wisdom, you can find a guide for just about
anything. Consider this your **ultra-condensed Docker cheat sheet** — the most
useful commands, tips, and gotchas that I’ve collected while wrangling
containers. Think of it as the espresso shot before the full Docker latte.

.. _docker-installation:

Installation
=============

* Ubuntu users: ``Docker`` is part of the latest Ubuntu versions 

* MacOs users: you can use the precompiled package directly from the `Docker webpage <https://www.docker.com/>`__.

Documentation
==================

For more details, the official source is always your friend:


* `Docker reference manual <https://docs.docker.com/>`_



Configuration
=============

At its core, a Docker container needs just two files to define its behavior:

* :file:`docker-compose.yml` — orchestrates services, networks, and volumes.

* :file:`Dockerfile` — describes how to build the container environment.




Useful tips
=============


Quick commands and tricks for keeping your Docker life sane.


Useful commands
---------------


* **List all volumes**:

    .. code:: bash

        docker volume ls

*Pro tip*: Volumes are where your data lives. Keep an eye on them.


How to change the default data directory
------------------------------------------

By default, Docker stores containers, images, and volumes in
:file:`/var/lib/docker`. If your root filesystem is limited, this may fill up
fast. Here’s how to move it to a new location (Ubuntu/Debian example):

(Extracted from `How to Change Docker’s Default Data Directory <https://linuxiac.com/how-to-change-docker-data-directory/>`_)

.. important:

    You will need sudo privileges for these steps.

1. Stop Docker:

    .. code:: bash

        sudo systemctl stop docker.service
        sudo systemctl stop docker.socket

2. Create a new location and move the data:

    .. code:: bash

        sudo mkdir /data/docker
        sudo mv /var/lib/docker /data/docker

3. Update  Docker configuration

    Edit the daemon config:

    .. code:: bash

        sudo vi /etc/docker/daemon.json 

    Add the following:

    .. code:: text

        { 
            "data-root": "/data/docker"
        }

4. Restart Docker:

    .. code:: bash

        sudo systemctl start docker.socket
        sudo systemctl start docker.service

5. Verify that Docker is running:


    .. code:: bash

        sudo systemctl status docker

    should output something like:

    .. code:: text

        ● docker.service - Docker Application Container Engine
        Loaded: loaded (/usr/lib/systemd/system/docker.service; enabled; preset: enabled)
        Active: active (running) since Sat 2025-10-18 11:05:22 EEST; 3h 24min ago
        Invocation: f2052aa83cfe491ebed5264d230aa966
        TriggeredBy: ● docker.socket
        Docs: https://docs.docker.com
        Main PID: 2353686 (dockerd)
        Tasks: 16
        Memory: 37.9M (peak: 41M)
        CPU: 38.638s
        CGroup: /system.slice/docker.service
             └─2353686 /usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock

        Oct 18 14:24:50 deer dockerd[2353686]: time="2025-10-18T14:24:50.662993087+03:00" level=info msg="No non-localhost DNS nameservers >
        Oct 18 14:24:50 deer dockerd[2353686]: time="2025-10-18T14:24:50.740593836+03:00" level=info msg="ignoring event" container=dae50ac>

6. Check that Docker is using the new directory:

    .. code:: bash

        $ sudo docker info | grep Root
        Docker Root Dir: /data/docker


Cleaning up space
----------------------

Docker loves to hoard. Let’s declutter.


(Extracted from `How to clear Docker cache and free up space on your system <https://depot.dev/blog/docker-clear-cache>`_)

* **Check disk usage**:

  .. code:: bash

      docker system df 

* **Remove stopped containers**:

  .. code:: bash

      docker container prune      # Removes all stopped containers

* **Remove all containers**:

  .. code:: bash

      docker stop $(docker ps -q) # Stop all containers
      docker container prune      # Remove stopped containers

* **Remove images**:

  .. code:: bash

      docker image prune -f       # Remove unused images
      docker image prune -a -f    # Remove all images

* **Remove volumes**:

  .. code:: bash

      docker volume prune -f      # Remove unused volumes
      docker volume rm -a         # Remove all volumes

* **Remove build cache**:

  .. code:: bash

      docker builder prune        # Old style
      docker buildx prune         # New style

* **Remove everything (use with caution!)**:

  .. code:: bash

      docker system prune -f

.. tip::

   Think of this like tidying your desk: removing containers and images you
   no longer need keeps Docker nimble, and saves your storage from quietly
   mutating into a monster.


