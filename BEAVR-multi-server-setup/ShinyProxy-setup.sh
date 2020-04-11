#!/bin/bash

#This script installs ShinyProxy on Ubuntu
#Script author: https://github.com/developerpiru/

echo ""
echo "=========================================================================="
echo "This script installs ShinyProxy on Ubuntu"
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
echo "Downloading ShinyProxy..."
wget https://www.shinyproxy.io/downloads/shinyproxy-2.3.0.jar

echo ""
echo "=========================================================================="
echo "ShinyProxy has been downloaded and Docker has been configured to use with" 
echo "ShinyProxy!"
echo ""
echo "Start a Docker container and then enter:"
echo ""
echo "java -jar shinyproxy-2.3.0.jar"
echo ""
echo "to start using your ShinyProxy server"
echo "=========================================================================="
