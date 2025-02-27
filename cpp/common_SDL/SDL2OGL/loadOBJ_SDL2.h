#include <iostream>
#include <fstream>
#include <string>

int loadObjToList( bool loadNormals, int maxn, char* filename ){
	using namespace std;
	Vec3f* vs  = new  Vec3f[maxn];
	Vec3f* vns;
	if( loadNormals ) vns = new Vec3f[4*maxn]; 
	int ilist = opengl1renderer.genLists(1);
	opengl1renderer.newList( ilist, GL_COMPILE );
	opengl1renderer.color3f( 0.9f, 0.9f, 0.9f );	
	//drawBox( -1, 1, -1, 1, -1, 1, 0.5, 0.5, 0.5 );
	int  x0 = -1,x1 = +1,y0 = -1,y1 = +1,z0 = -1,z1 = +1;
	ifstream myfile (filename);
	if ( myfile.is_open() )	{
		int istart,iend;
		string line;
		int iv  = 0;
		int ivn = 0;
		while ( getline (myfile,line) )		{
			// reading faces
			if( line.compare(0,2,"v ")==0 ){
				string word;
				istart   = line.find(" ",0     ); 
				iend     = line.find(" ",istart+1);
				word     = line.substr(istart,iend-istart);
				vs[iv].x = stof( word );
				istart   = iend;
				iend     = line.find(" ",istart+1);
				word     = line.substr(istart,iend-istart);
				vs[iv].y = stof( word );
				word     = line.substr(iend+1);
				vs[iv].z = stof( word );
				//printf( " _%f_%f_%f_ ",vs[3*iv],vs[3*iv+1],vs[3*iv+2] ); printf( line.c_str() ); printf( "\n" );
				iv++;
			} else
			// reading normals
			if( loadNormals && ( line.compare(0,2,"vn")==0 ) ){
				string word;
				istart    = line.find(" ",0     ); 
				iend      = line.find(" ",istart+1);
				word      = line.substr(istart,iend-istart);
				vns[ivn].x = stof( word );
				istart    = iend;
				iend      = line.find(" ",istart+1);
				word      = line.substr(istart,iend-istart);
				vns[ivn].y = stof( word );
				word      = line.substr(iend+1);
				vns[ivn].z = stof( word );
				//printf( " _%f_%f_%f_ ",vns[ivn].x,vns[ivn].y,vns[ivn].z ); printf( line.c_str() ); printf( "\n" );
				ivn++;
			} else
			// reading faces
			if( line.compare(0,2,"f ")==0 ){
				int fciv[4];
				int fcin[4];
				int nv=0;
				istart = line.find(" ",0 );
				iend   = line.find(" ",istart+1);
				while( iend != string::npos ){
					int islash  = line.find("/",istart+1);
					string word = line.substr(istart,islash-istart);
					fciv[nv]    = stoi( word );
					if (loadNormals){
						string word   = line.substr(islash+2,iend-islash-1);
						fcin[nv]      = stoi( word );
						//cout << word << "\n";
					}
					nv++;
					istart = iend;
					iend = line.find(" ",istart+1);
				}
				//cout << line << "\n";
				if (loadNormals){
				//if (false){
					opengl1renderer.enable (GL_LIGHTING);
					if      ( nv == 3 ){ opengl1renderer.begin(GL_TRIANGLES); }
					else if ( nv == 4 ){ opengl1renderer.begin(GL_QUADS);     };
					for ( int j=0; j<nv; j++ ){ 
						int i  = fciv[j]-1;
						int in = fcin[j]-1;
						opengl1renderer.normal3f(vns[in].x,vns[in].y,vns[in].z); 
						opengl1renderer.vertex3f( vs[i ].x, vs[i ].y, vs[i ].z );   
					}
					opengl1renderer.end();
/*
					opengl1renderer.disable (GL_LIGHTING);
					for ( int j=0; j<nv; j++ ){ 
						int i  = fciv[j]-1;
						int in = fcin[j]-1;
						opengl1renderer.begin(GL_LINES);
							opengl1renderer.vertex3f( vs[i ].x, vs[i ].y, vs[i ].z );  
							opengl1renderer.vertex3f( vs[i ].x + vns[in].x, vs[i ].y + vns[in].y, vs[i ].z + vns[in].z );
						opengl1renderer.end();
					}
*/
				}else{ // compute normals

					int   i0 = fciv[0]-1; int   i1 = fciv[1]-1; int   i2 = fciv[2]-1;
					Vec3f a       = vs[i0]-vs[i1];
					Vec3f b       = vs[i2]-vs[i1];
					Vec3f normal  = a.cross( b );
					normal       *= -1/sqrt(normal.mag2()); 
					opengl1renderer.enable (GL_LIGHTING);
					if      ( nv == 3 ){ opengl1renderer.begin(GL_TRIANGLES); }
					else if ( nv == 4 ){ opengl1renderer.begin(GL_QUADS); };
					for ( int j=0; j<nv; j++ ){ 
						int i = fciv[j]-1; 
						opengl1renderer.normal3f(normal.x,normal.y,normal.z); 
						opengl1renderer.vertex3f( vs[i].x, vs[i].y, vs[i].z );   
					}
					opengl1renderer.end();
/*
					int i = fciv[1]-1;
					opengl1renderer.disable (GL_LIGHTING);
					opengl1renderer.begin(GL_LINES);
						Vec3f center = vs[i1] + (a + b)*0.5f;
						opengl1renderer.vertex3f( center.x           , center.y           , center.z            );
						opengl1renderer.vertex3f( center.x + normal.x, center.y + normal.y, center.z + normal.z );
					opengl1renderer.end();
*/
				}
				//printf( " %i %f %f %f %f \n", i, nx,ny,nz, nx*nx+ny*ny+nz*nz );
				//cout << nv <<"           " << line << "\n";
				//break;
			}
		}
	myfile.close();
	}else cout << "Unable to open file"; 
	delete [] vs;
	if( loadNormals ) delete [] vns; 
	//delete [] fcs;
	opengl1renderer.endList();
	return ilist;
}

