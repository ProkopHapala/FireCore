/// @file browser.h
/// @brief Core directory browsing — file listing and extension filtering.
///
/// Browser reads a directory and separates entries into fileNames and subDirNames,
/// filtering files by registered extensions (e.g. .xyz, .mol, .mol2). Provides
/// checkExtension() and readDir() methods. Used as base class by BrowserView
/// (thumbnail browser) and other file-browsing components. No SDL or OpenGL deps.

#ifndef  browser_sdl_h
#define  browser_sdl_h

#include <vector>
#include <string>
#include <algorithm>
#include <dirent.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>

using namespace std;

class Browser{ public:
	enum JOB       { JOB_NONE,   JOB_THUMBS };
	enum MODE      { MODE_VIEW, MODE_THUMBS };
	enum TILETYPE  { TILE_NONE,  TILE_DIR, TILE_IMG };

	// files
	string work_dir;
    vector<string> fileNames;
    vector<string> subDirNames;
    unordered_map<string,int> extensions;

    //bool bMyDownScale = true;
    bool bMyDownScale = false;

	JOB  job;
	MODE mode;

	int curRow;
	int curCol;
	int curThumb;
	TILETYPE curType;

	Browser( string work_dir_ ){
		work_dir = work_dir_;
		fileNames    = vector<string>();
    	subDirNames = vector<string>();
	}

	bool checkExtension( const string& name ){
		int idot = name.find_last_of("."); 
		if(idot == (int)string::npos || idot == 0) return false;
		string ext = name.substr( idot + 1);
		//cout << name << " idot: " << idot << "  ext: " << ext << endl; 
        //printf( "ext: %s \n", ext.c_str() );
        auto found = extensions.find(ext);
		if( found != extensions.end() ){
            //printf( "found !!! \n" );
			return true;
		}
		return false;
	}

	int readDir ( const  string& dir ){
		fileNames.clear();
		subDirNames.clear();
		DIR *dp;
		struct dirent *entity;
		dp  = opendir(dir.c_str());
		if(  dp == NULL  ) {   
			printf( " Error %i opening %.100s\n", errno, dir.c_str() );  
			return -1;
		}else{
			subDirNames.push_back( string("..") );
			while ((entity = readdir(dp)) != NULL) {
				string fname = string(entity->d_name);
				if(fname[0]=='.' && (fname.size()==1 || (fname.size()==2 && fname[1]=='.'))) continue;
				bool isDir = false;
				bool isReg = false;
				string fpath = dir + "/" + fname;
				if(entity->d_type == DT_DIR) isDir = true;
				else if(entity->d_type == DT_REG) isReg = true;
				else {
					struct stat st;
					if(stat(fpath.c_str(), &st)==0){
						if(S_ISDIR(st.st_mode)) isDir = true;
						else if(S_ISREG(st.st_mode)) isReg = true;
					}
				}
				if(isDir){
                    if(fname[0]=='.') continue; // hide .vispy_mol_browser_cache etc.
                    subDirNames.push_back(fname);
                }
				else if(isReg && checkExtension(fname)) fileNames.push_back(fname);
			}
			closedir(dp);

			if(subDirNames.size() > 1)
				std::sort(subDirNames.begin()+1, subDirNames.end());
			std::sort(fileNames.begin(), fileNames.end());

			printf( "===== %i images in dir: %.100s \n", (int)fileNames.size(), dir.c_str() );
			for(int i=0; i<(int)fileNames.size(); i++){
				printf( " %i   %.100s\n", i, fileNames[i].c_str() ); 
			}
			printf( "===== %i dirs in dir: %.100s \n", (int)subDirNames.size(), dir.c_str() );
			for(int i=0; i<(int)subDirNames.size(); i++){
				printf( " %i   %.100s\n", i, subDirNames[i].c_str() ); 
			}

			return 1;
		}
	}

};

#endif