#!/bin/bash

#This script installs the OpenJDK Java environment on Ubuntu
#Script author: https://github.com/developerpiru/

echo ""
echo "=========================================================================="
echo "This script installs the OpenJDK Java environment on Ubuntu"
echo "Script author: https://github.com/developerpiru"
echo "=========================================================================="

echo ""
echo "Updating repositories..."
echo ""
sudo apt-get update -y

echo ""
echo "Installing Java OpenJDK..."
sudo apt-get install default-jdk -y

echo ""
echo "If you see the OpenJDK version number below, then OpenJDK has instlled successfully!"
echo ""
java -version
echo ""


