run:
	# To run on bash
	#docker run --rm -p 8787:8787 -v $(pwd)/src/:/home/rstudio/scripts rstudio
	docker run --rm -p 8787:8787 -v $(dir $(realpath $(firstword $(MAKEFILE_LIST))))src/:/home/rstudio/scripts rstudio
build:
	docker build -t rstudio .   
clean:
	docker images -f dangling=true
help:
	@echo start - Will start the docker container
	@echo build - Will build the docker image
	@echo clean - Will remove the hanging docker images