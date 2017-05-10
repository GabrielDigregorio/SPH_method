// Use: ./ppmToBath bathymetry.ppm bathymetry.dat
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
struct Ppm {
        int w;
        int h;
        int maxval;
        int* grid;
};
typedef struct Ppm Ppm;
int ppmReader(char* ppmName, Ppm* ppm, char* filename){
        // ---DEFINE ---//
        double xa=0.0;
        double xb=200.0; //beach: 200.0 //volcano: 700.0
        double ya=0.0;
        double yb=60.0; //beach: 60.0 //volcano: 700.0
        double max=30.0; //Maximum h value (white colour)
        //--------------//
        int i;
        double val;
        int gotWH=0;
        int gotMaxval=0;
        char buf[1000];
        char first[1], second[1];
        FILE *fp;
        FILE* fpnew=fopen(filename,"wb");
        fp = fopen(ppmName, "r");
        fgets(buf, sizeof(buf), fp);
        sscanf(buf, "%c%c", first, second);
        if(first[0]=='P' && (second [0]=='3'||second [0]=='5'||second [0]=='6')){ // Verifying file format
                while(gotWH==0 || gotMaxval==0){ // Reading w, h, maxval and skipping comments
                        if(gotWH==0){
                                fgets(buf, sizeof(buf), fp);
                                if(buf[0]!='#'){ // w and h value
                                        sscanf(buf, "%d %d", &(ppm->w), &(ppm->h));
                                        gotWH=1;
                                }
                        }
                        if(gotWH==1){
                                fgets(buf, sizeof(buf), fp);
                                if(buf[0]!='#'){ //maxval
                                        sscanf(buf, "%d", &(ppm->maxval));
                                        gotMaxval=1;
                                }
                        }
                }
                int w=ppm->w-1;
                int h=ppm->h-1;
                printf("\n %d %d \n", w,h);
                fwrite(&xa,sizeof(double),1, fpnew);
                fwrite(&xb,sizeof(double),1, fpnew);
                fwrite(&ya,sizeof(double),1, fpnew);
                fwrite(&yb,sizeof(double),1, fpnew);
                fwrite(&h,sizeof(int),1, fpnew);
                fwrite(&w,sizeof(int),1, fpnew);
                // Verifications
                assert((ppm->maxval)<=255 && "Wrong Maxval value.");
                // Allocating memory and reading the grid
                ppm->grid = (int*)malloc(3*(ppm->w)*(ppm->h)*sizeof(int));
                for (i=0; i<(3*(ppm->w)*(ppm->h)); ++i){
                        (ppm->grid)[i] = fgetc(fp);
                        if((i+1)%3==0){
                                val=((double)(255-(ppm->grid)[i]))/255.0*max;
                                fwrite(&val,sizeof(double),1, fpnew);
                        }
                }
                fclose(fp);
                fclose(fpnew);
                return 0;
        }
        else{
                printf("File is not in the .ppm format.\n");
                fclose(fp);
                return 1;
        }
}
void clear(Ppm* ppm){
        // Freeing the allocated grids
        free(ppm->grid);
}
int main(int argc, char **argv){
        if(argc==3){ // Single PPM
                Ppm ppm;
                printf("Processing %s ...\n", argv[1]);
                if(ppmReader(argv[1], &ppm, argv[2])==0){
                        clear(&ppm);
                }
        }
                        return 0;
}
