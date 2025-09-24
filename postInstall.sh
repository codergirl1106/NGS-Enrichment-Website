#!/bin/bash

# Fix any partially installed packages
dpkg --configure -a
apt-get install -f -y

# Hold conflicting ODBC packages so they aren't installed
apt-mark hold libodbc2 libodbcinst2 unixodbc-common

# Install minimal ODBC headers (won't trigger conflicts)
apt-get install -y unixodbc unixodbc-dev
