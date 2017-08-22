Instructions
============

To amend these images make changes in the appropriate docker file and build:

  docker build . -f <filename>

The unique sha associated with the id can be found from:

  docker images

This image should then be tagged appropriately, for example the images pushed to dockerhub are for example, pyneci-python27-ubuntu1704 :

  docker tag <sha id> "pynedist/meaningful_name" (the pynedist is the Dockerhub username)

If you are happy with the image, publish it to Dockerhub, first login as pynedist, then push

  docker login pynedist
  docker push pynedist/meaningful_name

If you have ammended image names, added python verison, there may need to be changes to the .travis.yml file.