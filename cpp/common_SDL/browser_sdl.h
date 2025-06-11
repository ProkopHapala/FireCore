#ifndef  browser_sdl_h
#define  browser_sdl_h


#include <SDL2/SDL.h>
#include "imageOps.h"


class BrowserSDL{
	public:
	enum JOB       { JOB_NONE,   JOB_THUMBS };
	enum MODE      { MODE_VIEW, MODE_THUMBS };
	enum TILETYPE  { TILE_NONE,  TILE_DIR, TILE_IMG };

	// files
	string work_dir;
    vector<string> imgNames;
    vector<string> subDirNames;
	// images
	int bd,tileW,tileH;
	SDL_Rect thumbRect;
	SDL_Rect thisRect;
	vector<SDL_Surface*>  thumbs;
	SDL_Surface         * thisImage;

    //bool bMyDownScale = true;
    bool bMyDownScale = false;


	JOB job;
	MODE mode;

	int iThumbCurr;
	int row0;

    int width;
    int height;

	int curRow;
	int curCol;
	int curThumb;
	TILETYPE curType;

	int nRowScreen;
	int nColScreen;

	int nRowImg;
	int nRowDir;

	Uint32 clrBg;
	Uint32 clrThumb;
	Uint32 clrDir;
	Uint32 clrCur;
	Uint32 clrText;

	Screen2d * screen;
	TTF_Font * font;

	Browser( string work_dir_ ){
		// files
		work_dir = work_dir_;
		imgNames    = vector<string>();
    	subDirNames = vector<string>();
		// images
		thumbs = vector<SDL_Surface*>();
	}

	void init(){
		clrBg     = 0xFF808080;
		clrThumb  = 0xFF606060;
		clrDir    = 0xFFF08000;
		clrCur    = 0xFF80FF00;
		clrText   = 0xFF00008F;
		thisImage = NULL;
		bd        = 2;
		tileW     = 180;
		tileH     = 180;
		thisRect.x   = 0;
		thisRect.y   = 0;
		thumbRect.x  = 0;
		thumbRect.y  = 0;
		thumbRect.w  = tileW-bd;
		thumbRect.h  = tileW-bd;
		nColScreen = screen->surface->w / tileW;
		nRowScreen = screen->surface->h / tileH;
		mode = MODE_THUMBS;
		job  = JOB_NONE;
		
	}

	bool isImage( const string& name ){
		int idot = name.find_last_of("."); 
		string ext = name.substr( idot + 1);
		//cout << name << " idot: " << idot << "  ext: " << ext << endl; 
		if( 
            (ext == "jpg") || (ext == "jpeg") || (ext == "JPG") || (ext == "JPEG") ||
            (ext == "png") || (ext == "PNG") || (ext == "bmp") || (ext == "BMP") || (ext == "tif") || (ext == "TIF") 
        ){
			return true;
		}
		return false;
	}

	int readDir ( const  string& dir ){
		imgNames.clear();
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
					if( isImage( fname ) ) {
						//cout << fname << endl;
						imgNames.push_back( fname );
					}
				}
				//std::cout << "Not a file or directory: " << entity->d_name << std::endl;
			}
			nRowImg = ceil(    imgNames.size() / float( nColScreen ) );
			nRowDir = ceil( subDirNames.size() / float( nColScreen ) );
			closedir(dp);

			printf( "===== %i images in dir: %.100s \n", (int)imgNames.size(), dir.c_str() );
			for(int i=0; i<(int)imgNames.size(); i++){
				printf( " %i   %.100s\n", i, imgNames[i].c_str() ); 
			}
			printf( "===== %i dirs in dir: %.100s \n", (int)subDirNames.size(), dir.c_str() );
			for(int i=0; i<(int)subDirNames.size(); i++){
				printf( " %i   %.100s\n", i, subDirNames[i].c_str() ); 
			}

			return 1;
		}
	}

/*
	int setAndReadDir( const  string& dir ){
		string dir_; 
		if( dir != ".." ){
			dir_ = work_dir+"/"+dir;
		}else{
			dir_ = dir;
		}
		printf( " new dir %s \n", dir_.c_str() );
		int err = chdir( dir_.c_str() );
		if( err != 0 ) {
			printf( " cannot change directory to %s \n", dir_.c_str() );
			return -1;
		}else{
			char  wd[512];
			char* wd_ = getcwd( wd, sizeof(wd) );
			printf( " new dir %s \n", wd_ );
			work_dir = string(wd_);
			return readDir( work_dir ); ;
		}
	}
*/

	int setAndReadDir( const  string& dir ){
		printf( " change to dir: %.100s \n", dir.c_str() );
		int err = chdir( dir.c_str() );
		if( err != 0 ) {
			printf( " cannot change directory to: %.100s \n", dir.c_str() );
			return -1;
		}else{
			char  wd[512];
			char* wd_ = getcwd( wd, sizeof(wd) );
			work_dir = string(wd_);
			printf( " work_dir: %.100s \n", work_dir.c_str() );
            SDL_SetWindowTitle( screen->window, work_dir.c_str() );
			return readDir( work_dir ); ;
		}
	}

	void deleteThumbs(){
		 while ( !thumbs.empty() ){
			SDL_FreeSurface( thumbs.back() );
			thumbs.pop_back();
		}
	}

	void readThumb( int iThumb ){
		string * fname = &imgNames[iThumb];
		thisImage = IMG_Load( fname->c_str() );
		if ( thisImage !=NULL ){
			//cout << " loaded: " << *fname << endl;
			printf( "readThumb %i loaded: %.100s \n", iThumb, fname->c_str() );
			SDL_Surface * thumb = SDL_CreateRGBSurface(0, thumbRect.w, thumbRect.h, 32, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );
			thisRect.w = thisImage->w;
			thisRect.h = thisImage->h;
			float aspect = thisRect.w/float(thisRect.h);
			SDL_Rect rdest;
			if( aspect < 1 ){
				rdest.h = thumbRect.h;
				rdest.w = int( thumbRect.w * aspect );
			}else{
				rdest.h = int( thumbRect.h / aspect );
				rdest.w = thumbRect.w;
			}
			rdest.x = (thumbRect.w - rdest.w) / 2;
			rdest.y = (thumbRect.h - rdest.h) / 2;
			
            if(bMyDownScale){
			    downScale( thisImage, thumb, &rdest );
            }else{
                SDL_BlitScaled( thisImage, &thisRect,  thumb, &rdest );
            }

			SDL_Color clr_txt;
 			clr_txt.a = ( clrText>>24 )&0xFF;
			clr_txt.r = ( clrText>>16 )&0xFF;
			clr_txt.g = ( clrText>>8  )&0xFF;
			clr_txt.b = ( clrText     )&0xFF;

            char txDebug[256];
            sprintf( txDebug, "nb %i", thisImage->pitch/thisImage->w );
			SDL_Surface* textSurf = TTF_RenderText_Solid( font, txDebug , clr_txt );
            //SDL_Surface* textSurf = TTF_RenderText_Solid( font, fname->c_str() , clr_txt );
            //modifyChanel(  thumb, 0 );

			SDL_BlitSurface( textSurf, &thumbRect, thumb, NULL );
			//SDL_BlitSurface( textSurf, &thumbRect, screen->surface, NULL );
			SDL_FreeSurface( textSurf );
			SDL_FreeSurface( thisImage );
			thumbs.push_back( thumb );
		} else{
			//cout << "cannot load: " << *fname << endl;
			printf( "readThumb %i cannot load: %.100s \n", iThumb, fname->c_str() );
			thumbs.push_back( NULL );
		}
	}

	void resetCur(){
		curRow   = 0;
		curCol   = 0;
		curThumb = 0;
		curType  = TILE_DIR;
	}

	int thumb2Col   ( int iThumb ){ return iThumb % nColScreen;	 }
	int thumb2Row   ( int iThumb ){ return iThumb / nColScreen;  }

	int rowCol2thumb( int row, int col ){ return row*nColScreen + col; }
	int rowCol2thumb( int row, int col, TILETYPE& tileType ){
		tileType = TILE_DIR;
		if( row >= nRowDir ){
			row -= nRowDir;
			tileType = TILE_IMG;
		}
		return rowCol2thumb( row, col ); 
	}

	void view( const string& fname ){
		//cout << " load: " << fname << endl;

		printf( "view loading: %.100s \n", fname.c_str() );
		thisImage = IMG_Load( fname.c_str() );
		if ( thisImage !=NULL ){
			//cout << " loaded: " << fname << endl;

            screen->resize( thisImage->w, thisImage->h );
            //SDL_SetWindowSize( screen->window, thisImage->w, thisImage->h ); 
            //screen->force_update();
            //surface = SDL_GetWindowSurface( window );
            //opengl1renderer.viewport(0, 0, windowWidth, windowHeight );

			printf( "loaded: %.100s \n", fname.c_str() );
			SDL_FillRect( screen->surface, NULL, clrBg );
			thisRect.w = thisImage->w;
			thisRect.h = thisImage->h;
			int sw = screen->surface->w;
			int sh = screen->surface->h;
			float aspect  = thisRect.w/float(thisRect.h);
			float saspect = sw/float(sh);
			SDL_Rect rdest;
			if( aspect < saspect ){
				rdest.h = sh;
				rdest.w = int( sw * aspect  / saspect  );
			}else{
				rdest.h = int( sh * saspect / aspect   );
				rdest.w = sw;
			}
			rdest.x = (screen->surface->w - rdest.w) / 2;
			rdest.y = (screen->surface->h - rdest.h) / 2;
			//SDL_BlitScaled( thisImage, &thisRect,  screen->surface, &rdest );
            if(bMyDownScale){
			    downScale( thisImage, screen->surface, &rdest );
            }else{
                SDL_BlitScaled( thisImage, &thisRect,  screen->surface, &rdest );
            }
		} else { 
			//cout << "cannot load: " << *fname << endl;
			printf( "view cannot load: %.100s \n", fname.c_str() );
			thumbs.push_back( NULL );
		}
	}

	void drawThumb( int iThumb, SDL_Surface * dest, Uint32 color ){
		//printf( " iThumb %i \n", iThumb );
		if ( iThumb < (int)thumbs.size() ){
			SDL_Surface * thumb = thumbs[iThumb];
			//printf( " iThumb %i \n", iThumb );
			if ( thumb !=NULL ){
				int irow = thumb2Row( iThumb ) + nRowDir - row0;
				if( (irow >=0) && (irow <= nRowScreen ) ){
					int icol = thumb2Col( iThumb ); 
					SDL_Rect rdest;
					rdest.x = icol*tileW;
					rdest.y = irow*tileH;
					rdest.w = thumbRect.w;
					rdest.h = thumbRect.h;
					SDL_FillRect( dest, &rdest, color );
					//printf( "  %i    %i %i \n", iThumb, thumb->w, thumb->h );
					SDL_BlitSurface( thumb, &thumbRect, dest, &rdest );
				}
			}
		}
	}

	void drawDir( int iThumb, SDL_Surface * dest, Uint32 color ){
		if ( iThumb < (int)subDirNames.size() ){
			string * fname = &subDirNames[iThumb];
			int irow = thumb2Row( iThumb ) - row0;
			if( (irow >=0) && (irow <= nRowScreen ) ){
				int icol = thumb2Col( iThumb ); 
				SDL_Rect rdest;
				rdest.x = icol*tileW;
				rdest.y = irow*tileH;
				rdest.w = thumbRect.w;
				rdest.h = thumbRect.h;
				SDL_FillRect( dest, &rdest, color );
				SDL_Color clr_txt;
	 			clr_txt.a = ( clrText>>24 )&0xFF;
				clr_txt.r = ( clrText>>16 )&0xFF;
				clr_txt.g = ( clrText>>8  )&0xFF;
				clr_txt.b = ( clrText     )&0xFF;
				SDL_Surface* textSurf = TTF_RenderText_Solid( font, fname->c_str() , clr_txt );
				SDL_BlitSurface( textSurf, &thumbRect, dest, &rdest );
				SDL_FreeSurface( textSurf );
			}
		}
	}


	void drawThumb( int iThumb ){
		int clr = clrThumb;
		if( iThumb == curThumb ) clr=clrCur;
		drawThumb( iThumb, screen->surface, clr );
	}

	void drawDir( int iThumb ){
		int clr = clrDir;
		if( iThumb == curThumb ) clr=clrCur;
		drawDir  ( iThumb, screen->surface, clr );
	}

	void drawTile( int iThumb, TILETYPE typ  ){
		if     ( typ == TILE_DIR ){ drawDir  (iThumb); }
		else if( typ == TILE_IMG ){ drawThumb(iThumb); }
	}

	void drawDirs( ){
		for( int iThumb=0; iThumb<(int)subDirNames.size(); iThumb++ ){
			drawDir( iThumb, screen->surface, clrDir );
		}
	}

	void drawThumbs( ){
		for( int iThumb=0; iThumb<(int)thumbs.size(); iThumb++ ){
			drawThumb( iThumb, screen->surface, clrThumb );
		}
	}

	void drawTiles( ){
		SDL_FillRect( screen->surface, NULL, clrBg );
		drawDirs( );
		drawThumbs( );
	}


	void setCol( int col ){
		if( ( col >= 0) && ( col < nColScreen ) ){
			TILETYPE oType = curType;
			int oThumb     = curThumb;
			curCol     = col;
			curThumb = rowCol2thumb( curRow, curCol, curType );
			drawTile(   oThumb, oType   );
			drawTile( curThumb, curType );
			//screen->updated = false;
			screen->force_update();
			printf( " row %i col %i  thumb %i  type %i \n ",  curRow, curCol, curThumb, curType );
		}
	}

	void setRow( int row ){
		if( ( row>= 0) && ( row<( nRowDir + nRowImg ) ) ){
			TILETYPE oType = curType;
			int oThumb     = curThumb;
			curRow = row;
			curThumb = rowCol2thumb( curRow, curCol, curType );
			drawTile(  oThumb, oType  );
			int srow =  curRow - row0;
			if  ( srow < 0 ){
				row0 = curRow;
				drawTiles();
			} else if( srow > nRowScreen ){
				row0 = curRow - nRowScreen;
				drawTiles();
			} else{
				drawTile( curThumb, curType );
			}
			//screen->updated = false;
			screen->force_update();
			printf( " row %i col %i  thumb %i  type %i \n ",  curRow, curCol, curThumb, curType );
		}
		
	}


	void moveCol( int d ){	
		if( mode == MODE_THUMBS ) setCol( curCol + d );	
		if( mode == MODE_VIEW ) {
			int newThumb = curThumb + d;
			if( (newThumb>0)&&( newThumb<(int)thumbs.size() ) ){
				curThumb = newThumb;
				SDL_FreeSurface( thisImage );
				view( imgNames[curThumb] );
				//screen->updated = false;
				screen->force_update();
			}
		};	
	}

	void moveRow( int d ){	
		if( mode == MODE_THUMBS ) setRow( curRow + d );	
	}

	void enter(){
		printf( "return \n" );
		if( mode == MODE_VIEW ){
            //SDL_SetWindowSize( screen->window, width, height );
            //screen->force_update();
            screen->resize( width, height );
			mode = MODE_THUMBS;
			SDL_FreeSurface( thisImage );
			drawTiles();
			drawTile( curThumb, curType );
			//screen->updated = false;
			screen->force_update();
		}else if ( mode == MODE_THUMBS ) {
			if       ( curType == TILE_IMG ){
				if ( curThumb < (int)imgNames.size()  ){
					mode = MODE_VIEW;
					view( imgNames[curThumb] );
				}
			}else if ( curType == TILE_DIR ){
				//if ( setAndReadDir( work_dir+"/"+subDirNames[curThumb] ) > 0){
				if ( curThumb < (int)subDirNames.size()  ){
					if ( setAndReadDir( subDirNames[curThumb] ) > 0){
						setJob_thumbs( );
					}
				}
			}
			//screen->updated = false;
			screen->force_update();
		}
	}

    void leave(){
		printf( "return \n" );
		if( mode == MODE_VIEW ){
		}else if ( mode == MODE_THUMBS ) {
            //if ( setAndReadDir( subDirNames[curThumb] ) > 0){
            if ( setAndReadDir( ".." ) > 0){
                setJob_thumbs( );
            }
			screen->force_update();
		}
	}

	void setJob_thumbs( ){
		iThumbCurr=0;
		job = JOB_THUMBS;
		deleteThumbs();
		SDL_FillRect( screen->surface, NULL, clrBg );
		drawDirs( );
		screen->updated = false;
	}

	void jobStep_thumbs(){
		if( iThumbCurr >= (int)imgNames.size() ){ 
			resetCur();
			drawTile( curThumb, curType );
			job=JOB_NONE;
			return;
		}else{
			readThumb( iThumbCurr );
			drawThumb( iThumbCurr, screen->surface, clrThumb );
			screen->updated = false;
			iThumbCurr++;
		}
	}

	void update(){
		switch( job ){
			case JOB_THUMBS:
				jobStep_thumbs();
				break;
			case JOB_NONE:
				break;
		}		
	}

	void close(){
		deleteThumbs();
		if( thisImage != NULL ) SDL_FreeSurface( thisImage );
		thisImage = NULL;
	}

};

#endif