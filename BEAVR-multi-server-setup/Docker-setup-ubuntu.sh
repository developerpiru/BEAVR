#!/bin/bash

# This script will install the latest verson of Docker Engine CE on Ubuntu
# Script author: https://github.com/developerpiru

echo ""
echo "=========================================================================="
echo "This script will install the latest verson of Docker Engine CE on Ubuntu"
echo "Script author: https://github.com/developerpiru"
echo "=========================================================================="

echo ""
echo "Updating repositories..."
echo ""
sudo apt-get update -y

echo ""
echo "Removing previous versions of Docker..."
echo ""
sudo apt-get remove docker docker-engine docker.io -y

echo ""
echo "Installing Docker..."
echo ""
sudo apt install docker.io -y

echo ""
echo "Configuring Docker to start on boot..."
echo ""
sudo systemctl start docker
sudo systemctl enable docker

echo ""
echo "If you see the Docker version number below, then Docker has instlled successfully!"
echo ""
docker --version
echo ""

