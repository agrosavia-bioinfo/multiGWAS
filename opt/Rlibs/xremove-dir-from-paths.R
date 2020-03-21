#!/usr/bin/Rscriot

args = commandArgs(trailingOnly = TRUE)

dirs = list.dirs (".", recursive=F)
currentDir = getwd()

for (dir in dirs) {
	libDirs = list.dirs (dir, recursive=F, full.names=F)
	if ("doc" %in% libDirs) {
		cmm = sprintf ("rm -rf %s/%s", dir, "doc")
		message (cmm)
		system (cmm)
	}
}


