#!/bin/bash

ulimit -n 4096
exec "$@"
