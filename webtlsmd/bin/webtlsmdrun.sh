#!/bin/bash

export TLSMD_ROOT=/home/tlsmd/tlsmd
export LOG_ROOT=$TLSMD_ROOT/log

## start web-app database server daemon
$TLSMD_ROOT/bin/webtlsmdd.py     > $LOG_ROOT/webtlsmdd.log    2>&1 &
sleep 1

## start the web-app tlsmd job running daemon
$TLSMD_ROOT/bin/webtlsmdrund.py  > $LOG_ROOT/webtlsmdrund.log 2>&1 &
