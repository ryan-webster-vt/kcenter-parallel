FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    openmpi-bin \
    libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Set work directory
WORKDIR /app

# Copy your code into container
COPY . /app

# Build your program (change filename if needed)
RUN mpicc mpi_kcenter.c -o kcenter -lm

# Default command (run with 4 processes as example)
CMD ["mpirun", "-np", "4", "./kcenter", "data.txt", "3", "1000", "42"]
