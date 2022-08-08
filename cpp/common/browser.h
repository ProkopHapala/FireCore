#ifndef  browser_sdl_h
#define  browser_sdl_h

#include <vector>
#include <string>
#include <dirent.h>
#include <stdio.h>
#include <unistd.h>

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
			//cout << "Error(" << errno << ") opening " << dir << endl; 
			printf( " Error %i opening %.100s\n", errno, dir.c_str() );  
			return -1;
		}else{
			subDirNames.push_back( string("..") );
			while ((entity = readdir(dp)) != NULL) {        
			   if(entity->d_type == DT_DIR){
					if(entity->d_name[0] != '.'){ // ignore '.' and '..'
						subDirNames.push_back( string(entity->d_name) );
					}
				}
				string fname = string(entity->d_name);
				if(entity->d_type == DT_REG){
					if( checkExtension( fname ) ) {
						//cout << fname << endl;
						fileNames.push_back( fname );
					}
				}
				//std::cout << "Not a file or directory: " << entity->d_name << std::endl;
			}
			//nRowImg = ceil(    fileNames.size() / float( nColScreen ) );
			//nRowDir = ceil( subDirNames.size() / float( nColScreen ) );
			closedir(dp);

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