#!/bin/bash

#This script installs ShinyProxy on Ubuntu and configures it for BEAVR
#Script author: https://github.com/developerpiru/

echo ""
echo "=========================================================================="
echo "This script installs ShinyProxy on Ubuntu and configures it for BEAVR"
echo "Script author: https://github.com/developerpiru"
echo "=========================================================================="

echo ""
echo "Configuring Docker proxy..."
sudo mkdir -p /etc/systemd/system/docker.service.d
sudo cp override.conf /etc/systemd/system/docker.service.d/override.conf

echo ""
echo "Restarting Docker..."
sudo systemctl daemon-reload
sudo systemctl restart docker

echo ""
echo "Downloading latest BEAVR container from Docker Hub..."
sudo docker pull pirunthan/beavr:latest

echo ""
echo "Downloading ShinyProxy..."
wget https://www.shinyproxy.io/downloads/shinyproxy-2.3.0.jar

echo ""
echo "=========================================================================="
echo "ShinyProxy and Docker have been configured to use with BEAVR"
echo ""
echo "Enter the following command to start the server:"
echo "java -jar shinyproxy-2.3.0.jar"
echo ""
echo "You can access the server from a browser at IP-ADDRESS:8080"
echo "See the application.yaml file for configuration options"
echo "=========================================================================="
