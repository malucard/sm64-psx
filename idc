#!/bin/sh

if [ -z "$MANAGER" ] || ! which $MANAGER; then
	if which docker 2>/dev/null; then
		MANAGER=docker
	elif which podman 2>/dev/null; then
		MANAGER=podman
	else
		>&2 echo "docker or podman not found"
		exit 1
	fi
fi

img=idc-tmp-$RANDOM
if $MANAGER build --rm -t $img .; then
	$MANAGER run -v .:/project -w /project --rm -it $img "$@"
	$MANAGER rmi $img
fi
