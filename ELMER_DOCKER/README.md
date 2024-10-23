# Project Description

This project uses Docker to containerize Elmer/Ice in an ubuntu environment. The Makefile provided allows you to easily build and run the Docker container.

## Prerequisites

- Docker: Ensure you have Docker installed on your machine. You can download and install Docker from [here](https://www.docker.com/products/docker-desktop).

## Installation

> ⚠️ **Warning:** Ensure Docker is running before executing the `make` commands.

1. **Install Docker**

   Follow the instructions on the Docker website to install Docker for your operating system.

2. **Open Terminal**

   Open your Terminal application and navigate to your the `ELMER_DOCKER` directory.

## Building the Docker Image

To build the Docker image, run the following command in your terminal:

```sh
make build
```

This command uses the `docker build` command to create a Docker image named `elmer`.

This command runs the Docker container and mounts a directory from your host machine to the container. By default, it mounts `/Users/cmosbeux` to `/home/elmer_usr/` inside the container.

> ⚠️ **Warning:** Change `/Users/cmosbeux` to the directory you want to access from the container. For example:

```sh
docker run -v /path/to/your/directory:/home/elmer_usr/ -it elmer
```

## Running the Docker Image

To run the Docker container, use the following command:

```sh
make run
```


