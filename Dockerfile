# Use the Cassiopee image on Docker Hub
FROM cassiopee486/cassiopee:main

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Set environment variables
ENV CASSIOPEE=/Cassiopee
ENV FAST=/Fast
ENV MACHINE=ubuntu

# Do not prevent mpirun to run as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Set the working directory in the container
WORKDIR $CASSIOPEE

# Copy the current directory contents into the container
COPY . $FAST

# Source environment and run install script
RUN . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 \
    && cd $FAST/Fast \
    && ./install

# Change the default shell to be the bash shell
SHELL ["/bin/bash", "-c"] 

# Define the default command to run the application: start an interactive shell
# session
ENTRYPOINT . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 && /bin/bash
