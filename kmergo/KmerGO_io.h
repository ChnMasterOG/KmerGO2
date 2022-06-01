/*
  Project Name	: KmerGO
  Version		: 2.0.0
  Author		: Qi Chen
  Date			: 2022-05-30
*/

#ifndef _KEMRGO_IO_H
#define _KMERGO_IO_H

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include <io.h>
#include <direct.h>
#elif defined(__linux__) || defined(linux)
#include <sys/io.h>
#include <dirent.h>
#include <unistd.h>
#endif

#include <iostream>
#include <vector>
#include <stdio.h>
#include <string.h>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#define popen	_popen
#define pclose	_pclose
#elif defined(__linux__) || defined(linux)
#define popen	popen
#define pclose	pclose
#endif

inline std::string joinPath(std::string path1, std::string path2)
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	/* replace / to \ */
	while (1)
	{
		std::string::size_type position;
		position = path1.find("/");
		if (position == path1.npos) break;
		path1 = path1.replace(position, 1, "\\");
	}
	while (1)
	{
		std::string::size_type position;
		position = path2.find("/");
		if (position == path2.npos) break;
		path2 = path2.replace(position, 1, "\\");
	}
	/* join two paths */
	std::string::iterator it = path1.end();
	--it;
	if (*it != '\\') path1.insert(path1.end(), '\\');
#elif defined(__linux__) || defined(linux)
	/* replace \ to / */
	while (1)
	{
		std::string::size_type position;
		position = path1.find("\\");
		if (position == path1.npos) break;
		path1 = path1.replace(position, 1, "/");
	}
	while (1)
	{
		std::string::size_type position;
		position = path2.find("\\");
		if (position == path2.npos) break;
		path2 = path2.replace(position, 1, "/");
	}
	/* join two paths */
	std::string::iterator it = path1.end();
	--it;
	if (*it != '/') path1.insert(path1.end(), '/');
#endif
	path1.insert(path1.size(), path2);
	return path1;
}

inline void getFiles(std::string path, std::vector<std::string>& files, bool absolute_path)
{
	char path_str[1024];
	strcpy(path_str, path.c_str());
	files.clear();

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	char path_buf[1024] = "";	// absolute path

	/* get real path */
	_fullpath(path_buf, path_str, 1024);
	if (path_buf[strlen(path_buf) - 1] != '\\')
	{
		path_buf[strlen(path_buf)] = '\\';
		path_buf[strlen(path_buf) + 1] = '\0';
	}
	/* get files */
	std::string p;
	int64 hFile = 0;
	struct _finddata_t fileinfo;
	if ((hFile = _findfirst(p.assign(std::string(path_buf)).append("*").c_str(), &fileinfo)) != -1) {
		do
		{
			if (fileinfo.name[0] == '.') continue;
			if (absolute_path == true) files.push_back(p.assign(std::string(path_buf)).append(fileinfo.name));
			else files.push_back(fileinfo.name);
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
#elif defined(__linux__) || defined(linux)
	/* get files */
	DIR* dp;
	struct dirent* dirp;
	if ((dp = opendir(path_str)) != NULL)
	{
		while ((dirp = readdir(dp)) != NULL) {
			if (dirp->d_name[0] == '.') continue;
			if (absolute_path == true) files.push_back(joinPath(path_str, dirp->d_name));
			else files.push_back(dirp->d_name);
		}
	}
#endif
}

inline std::string getCWD()
{
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	return std::string(_getcwd(NULL, 0));
#elif defined(__linux__) || defined(linux)
	return std::string(getcwd(NULL, 0));
#endif
}

inline uint64 getFileSize(std::string fpath)
{
	std::ifstream in(fpath, std::ifstream::ate | std::ifstream::binary);
	return in.tellg();
}

#endif
