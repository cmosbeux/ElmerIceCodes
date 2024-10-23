# Makefile

DOCKER_IMAGE_NAME := elmer
DOCKER_CONTAINER_NAME := elmer-container

build:
	#docker build --build-arg SOURCE_DIR=$(ELMER_SOURCE) -t $(DOCKER_IMAGE_NAME) .
	docker build -t $(DOCKER_IMAGE_NAME) .

run:
	docker run -v /Users/cmosbeux:/home/elmer_usr/ -it $(DOCKER_IMAGE_NAME)

link:
	@if [ -z "$(DIR)" ] ; then \
		rm -f elmer ; \
	else \
		ln -sfn "$(DIR)" elmer ; \
	fi

.PHONY: build run link

