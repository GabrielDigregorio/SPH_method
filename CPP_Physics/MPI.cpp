///**************************************************************************
/// SOURCE: Functions related to MPI.
///**************************************************************************
#include "Main.h"
#include "Interface.h"
#include "Physics.h"
#include "Tools.h"
#include "Structures.h"

// For particle exchange/sharing
enum overlap {leftOverlap, noOverlap, rightOverlap, NB_OVERLAP_VALUE};
enum migrate {noMigrate, leftMigrate, rightMigrate, NB_MIGRATE_VALUE};
enum mpiMessage {sendOverlap, recvOverlap, sendMigrate, recvMigrate, sendData, recvData, NB_MPIMESSAGE_VALUE};


/*
Input:
    - globalField: filled only on process #0, to be scattered
    - localField: partial field specific to each process. Will contain
    both the domain to be solved and the halos. Need to be done in two
    steps because MPI_Scatterv does not allow overlaps.
    - parameter: to get kh and integrationMethod
Ouput:
    - errorFlag: tells if the number of processor was acceptable or not
*/
void scatterField(Field* globalField, Field* localField, Parameter* parameter,
    SubdomainInfo &subdomainInfo){
    // Basic MPI process information
    int nTasks = subdomainInfo.nTasks;
    int procID = subdomainInfo.procID;
    // Box size (bigger if RK2 to avoid sorting twice at each time step)
    double boxSize = boxSizeCalc(parameter->kh, parameter->integrationMethod);
    // Checks if the number of processor appropriate (only node 0 has the info)
    int nTotalBoxesX;
    if(procID == 0){
        for(int i=0 ; i<3 ; i++){
            localField->l[i] = globalField->l[i]; // x component will be changed
            localField->u[i] = globalField->u[i]; // x component will be changed
        }
        nTotalBoxesX = ceil((globalField->u[0] - globalField->l[0])/boxSize);
        if(nTotalBoxesX < 3*nTasks){
            std::cout << "Too much processors for the domain" << std::endl;
            std::cout << "-> not yet handled !!!" << std::endl;
        }
        else{std::cout << "Appropriate number of processors" << std::endl;}
    }
    // Broadcasts the total number of boxes along x
    MPI_Bcast(&nTotalBoxesX, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> startBoxX(nTasks+1); // The last element helps the last process
    for(int i=0 ; i <= nTasks ; i++){
        startBoxX[i] = (nTotalBoxesX*i) / nTasks; // TO OPTIMIZE ???
    }

    // Broadcasts the global l and u and determines the correct l[0] and u[0]
    MPI_Bcast(localField->l, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(localField->u, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double globall0 = localField->l[0];

    // Number of boxes along the Y and Z directions
    int nBoxesY = ceil((localField->u[1] - localField->l[1])/boxSize);
    int nBoxesZ = ceil((localField->u[2] - localField->l[2])/boxSize);

    // Left boundary and starting box
    if(procID == 0){
        localField->l[0] = globall0;
        subdomainInfo.startingBox = 0;
    }
    else{
        localField->l[0] = globall0 + (startBoxX[procID]-1)*boxSize;
        subdomainInfo.startingBox = nBoxesY*nBoxesZ;
    }
    // Right boundary and ending box
    if(procID == nTasks - 1){
        localField->u[0] = globall0 + startBoxX[procID+1]*boxSize;
        subdomainInfo.endingBox = subdomainInfo.startingBox + (startBoxX[procID+1]-startBoxX[procID])*nBoxesY*nBoxesZ;
    }
    else{
        localField->u[0] = globall0 + (startBoxX[procID+1]+1) *boxSize;
        subdomainInfo.endingBox = subdomainInfo.startingBox + (startBoxX[procID+1]-startBoxX[procID])*nBoxesY*nBoxesZ;
    }

    std::cout << procID << ": l and u: " << localField->l[0] << " " << localField->u[0] << std::endl;
    std::cout << procID << ": boxes  : " << subdomainInfo.startingBox << " " << subdomainInfo.endingBox << std::endl;


    // Computes indices and sorts particles
    std::vector<int> nPartNode(nTasks, 0);
    std::vector< std::pair<int,int> > domainIndex;
    std::vector<double> limits(nTasks); // left boundary of each domain along x
    std::vector<int> offset(nTasks);

    if(procID==0){
        for(int i=0 ; i<nTasks ; i++)
            limits[i] = globall0 + startBoxX[i]*boxSize; // Without overlap
        computeDomainIndex(globalField->pos[0], limits, nPartNode, domainIndex, nTasks);
        sortParticles(*globalField, domainIndex); // BE SURE NOT COPY IS PERFORMED !!!!!!!!!!
        // Offset vector
        offset[0] = 0;
        for(int i=1 ; i<nTasks ; i++){
            offset[i] = offset[i-1]+nPartNode[i-1];
        }
        /* DISPLAY
        for(int i=0; i<domainIndex.size(); ++i){
            std::cout<<domainIndex[i].first <<" "<< domainIndex[i].second <<std::endl;
        }
        for(int i=0 ; i<nTasks ; i++)
            std::cout << nPartNode[i] << " ";
        std::cout << std::endl;
        //*/
    }


    // Shares the number of particle per domain and prepares the vector size
    MPI_Scatter(&nPartNode[0], 1, MPI_INT, &localField->nTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i=0 ; i<3 ; i++){
        localField->pos[i].resize(localField->nTotal);
        localField->speed[i].resize(localField->nTotal);
    }
    localField->density.resize(localField->nTotal);
    localField->pressure.resize(localField->nTotal);
    localField->mass.resize(localField->nTotal);
    localField->type.resize(localField->nTotal);

    std::cout << procID << ": "<<localField->nTotal << std::endl;

    // Scatters globalField into localFields
    for(int i=0 ; i<3 ; i++){
        MPI_Scatterv(&(globalField->pos[i][0]), &nPartNode[0], &offset[0], MPI_DOUBLE,
                &(localField->pos[i][0]), localField->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&globalField->speed[i][0], &nPartNode[0], &offset[0], MPI_DOUBLE,
                &localField->speed[i][0], localField->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    MPI_Scatterv(&(globalField->density[0]), &nPartNode[0], &offset[0], MPI_DOUBLE,
            &(localField->density[0]), localField->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&globalField->pressure[0], &nPartNode[0], &offset[0], MPI_DOUBLE,
            &localField->pressure[0], localField->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&globalField->mass[0], &nPartNode[0], &offset[0], MPI_DOUBLE,
            &localField->mass[0], localField->nTotal, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(globalField->type[0]), &nPartNode[0], &offset[0], MPI_INT,
            &(localField->type[0]), localField->nTotal, MPI_INT, 0, MPI_COMM_WORLD);

    // Computes nFree, nMoving and nFixed (IS IT REALLY USEFUL ??)
    localField->nFree = 0;
    localField->nFixed = 0;
    localField->nMoving = 0;
    for(int i=0 ; i<localField->nTotal ; i++){
        switch(localField->type[i]){
            case freePart:
            localField->nFree++;
            break;
            case fixedPart:
            localField->nFixed++;
            break;
            default:
            localField->nMoving++;
        }
    }

    std::cout << "Processor " << procID << " has " << localField->nFree << " free particles, " << localField->nFixed << " fixed particles and " << localField->nMoving << " moving particles.\n" << std::endl;

    // Sharing boundaries
    //shareOverlap(&localField, boxSize, procID, nTasks);

    // NEED TO DETERMINE STARTINGPARTICLE AND ENDINGPARTICLE
    // !!!!!

}

/* Gathers all the current fields into the global Field */
void gatherField(Field* globalField, Field* localField, SubdomainInfo &subdomainInfo){
    // Gathers the number of particles for each node in node 0
    int nbPart = subdomainInfo.endingParticle - subdomainInfo.startingParticle;
    std::vector<int> allNbPart;
    if(subdomainInfo.procID == 0)
        allNbPart.resize(subdomainInfo.nTasks);
    MPI_Gather(&nbPart, 1, MPI_INT, &(allNbPart[0]), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Computes the offsets for scatterv
    std::vector<int> offsets;
    if(subdomainInfo.procID == 0){
        offsets.resize(subdomainInfo.nTasks);
        offsets[0] = 0;
        for(int i = 1 ; i<subdomainInfo.nTasks ; i++)
            offsets[i] = offsets[i-1] + allNbPart[i-1];
    }

    // Gathers the fields
    for(int i=0 ; i<3 ; i++){
        MPI_Gatherv(&(localField->pos[i][0]), nbPart, MPI_DOUBLE,
                &(globalField->pos[i][0]), &(allNbPart[0]), &(offsets[0]), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
        MPI_Gatherv(&(localField->speed[i][0]), nbPart, MPI_DOUBLE,
                &(globalField->speed[i][0]), &(allNbPart[0]), &(offsets[0]), MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    }
    MPI_Gatherv(&(localField->density[0]), nbPart, MPI_DOUBLE,
            &(globalField->density[0]), &(allNbPart[0]), &(offsets[0]), MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    MPI_Gatherv(&(localField->pressure[0]), nbPart, MPI_DOUBLE,
            &(globalField->pressure[0]), &(allNbPart[0]), &(offsets[0]), MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    MPI_Gatherv(&(localField->mass[0]), nbPart, MPI_DOUBLE,
            &(globalField->mass[0]), &(allNbPart[0]), &(offsets[0]), MPI_DOUBLE,
            0, MPI_COMM_WORLD);
    MPI_Gatherv(&(localField->type[0]), nbPart, MPI_INT,
            &(globalField->type[0]), &(allNbPart[0]), &(offsets[0]), MPI_INT,
            0, MPI_COMM_WORLD);



    // INFORMATION
    if(subdomainInfo.procID == 0){
        std::cout << "\n\n\nGather information: " << std::endl;
        std::cout << "Number of particles per domain: " << std::endl;
        for(int i = 0 ; i<subdomainInfo.nTasks ; i++)
            std::cout << allNbPart[i] << " ";
        std::cout << "\nOffset for each domain: " << std::endl;
        for(int i = 0 ; i<subdomainInfo.nTasks ; i++)
            std::cout << offsets[i] << " ";
        std::cout << std::endl;
    }

}





void gatherField(Field* globalField, Field* localField){
    // Gather all the current fields into the global Field (for output file writing for example)

    globalField->nFree = localField->nFree;
    globalField->nFixed = localField->nFixed;
    globalField->nMoving = localField->nMoving;
    globalField->nTotal = localField->nTotal;
    for (int i = 0; i < 3; i++){
      globalField->l[i] = localField->l[i];
      globalField->u[i] = localField->u[i];
      globalField->pos[i] = localField->pos[i];
      globalField->speed[i] = localField->speed[i];
    }
    globalField->mass = localField->mass;
    globalField->type = localField->type;
    globalField->density = localField->density;
    globalField->pressure = localField->pressure;


}

void shareOverlap(Field& field, double u0, double l0, double boxSize, int procID, int nTasks){
    MPI_Status *status = NULL;
    std::vector< std::pair<int,int> >  indexOverlap;
    std::vector<double> recvVectL, recvVectR;
    int nOverlap[2];
    int sizeRecvOverlap;
    // Different possible cases (trivial if only one node)
    if(nTasks > 1){
        if(procID == 0){ // No left neighbor domain
            // inner domain: [l[0] , u[0] - 2*boxSize[
            // right edge: [u[0] - 2*boxSize , u[0] - boxSize[
            // right halo: [u[0] - boxSize , u[0]]
            computeOverlapIndex(field.pos[0],indexOverlap, nOverlap,
            l0+boxSize, l0+boxSize, u0-2*boxSize, u0-boxSize); // No Left Overlap
            sortParticles(field, indexOverlap);

            // Sends the edge to the right
            MPI_Send(&nOverlap[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
            MPI_Send(&field.pos[0], nOverlap[1], MPI_DOUBLE, procID+1, sendOverlap, MPI_COMM_WORLD);

            // Receives the halo from the right
            MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
            MPI_Recv(&recvVectR, sizeRecvOverlap, MPI_DOUBLE, procID+1, recvOverlap, MPI_COMM_WORLD, status);
            // Adding the particles to the field
            field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectL.end());
        }
        else if(procID == nTasks - 1){ // No right neighbor domain
            // left halo: [l[0] , l[0] + boxSize[
            // left edge: [l[0] + boxSize , l[0] + 2*boxSize[
            // inner domain: [l[0] + 2*boxSize , u[0][
            computeOverlapIndex(field.pos[0],indexOverlap, nOverlap,
            l0+boxSize, l0+2*boxSize, u0-boxSize, u0-boxSize); // No right Overlap
            sortParticles(field, indexOverlap);

            if(procID%2==0){
                // Sends the edge to the left (procID is even)
                MPI_Send(&nOverlap[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[0], MPI_DOUBLE, procID-1, sendOverlap, MPI_COMM_WORLD);

                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvOverlap, MPI_DOUBLE, procID-1, recvOverlap, MPI_COMM_WORLD, status);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].begin(), recvVectL.begin(), recvVectL.end());
           }else{
                // Receives the edge from the left (procID is odd)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvOverlap, MPI_DOUBLE, procID-1, recvOverlap, MPI_COMM_WORLD, status);

                // Sends the edge to the left (procID is odd)
                MPI_Send(&nOverlap[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[0], MPI_DOUBLE, procID-1, sendOverlap, MPI_COMM_WORLD);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].begin(), recvVectL.begin(), recvVectL.end());
            }
            // Receives/Send the halo from/to the left (procID is even/odd)

        }
        else{ // Domain in the middle
            // left halo: [l[0] , l[0] + boxSize[
            // left edge: [l[0] + boxSize , l[0] + 2*boxSize[
            // inner domain: [l[0] + 2*boxSize , u[0] - 2*boxSize[
            // right edge: [u[0] - 2*boxSize , u[0] - boxSize[
            // right halo: [u[0] - boxSize , u[0]]
            computeOverlapIndex(field.pos[0],indexOverlap, nOverlap,
            l0+boxSize, l0+2*boxSize, u0-2*boxSize, u0-boxSize);
            sortParticles(field, indexOverlap);
            if(procID%2==0){
                // Sends the edge to the left (procID is even)
                MPI_Send(&nOverlap[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[0], MPI_DOUBLE, procID-1, sendOverlap, MPI_COMM_WORLD);
                // Sends the edge to the right (procID is even)
                MPI_Send(&nOverlap[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[1], MPI_DOUBLE, procID+1, sendOverlap, MPI_COMM_WORLD);

                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvOverlap, MPI_DOUBLE, procID-1, recvOverlap, MPI_COMM_WORLD, status);
                // Receives the edge from the right (procID is even)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectR, sizeRecvOverlap, MPI_DOUBLE, procID+1, recvOverlap, MPI_COMM_WORLD, status);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].begin(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }else{
                // Receives the edge from the right (procID is even)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectR, sizeRecvOverlap, MPI_DOUBLE, procID+1, recvOverlap, MPI_COMM_WORLD, status);
                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvOverlap, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvOverlap, MPI_DOUBLE, procID-1, recvOverlap, MPI_COMM_WORLD, status);

                // Sends the edge to the right (procID is even)
                MPI_Send(&nOverlap[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[1], MPI_DOUBLE, procID+1, sendOverlap, MPI_COMM_WORLD);
                // Sends the edge to the left (procID is even)
                MPI_Send(&nOverlap[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nOverlap[0], MPI_DOUBLE, procID-1, sendOverlap, MPI_COMM_WORLD);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].begin(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }
        }
    }
    else{
        // Single processor
        return;
    }

}

void shareMigrate(Field& field, double u0, double l0, double boxSize, int procID, int nTasks){
    // Different possible cases (trivial if only one node)
    MPI_Status *status = NULL;
    std::vector< std::pair<int,int> >  indexMigrate;
    std::vector<double> recvVectL, recvVectR;
    int nMigrate[2];
    int sizeRecvMigrate;
    if(nTasks > 1){
        if(procID == 0){ // No left neighbor domain
            // inner domain: [l[0] , u[0] - 2*boxSize[
            // right edge: [u[0] - 2*boxSize , u[0] - boxSize[
            // right halo: [u[0] - boxSize , u[0]]
            computeMigrateIndex(field.pos[0], indexMigrate, nMigrate, l0, u0-boxSize);
            sortParticles(field, indexMigrate);

            // Sends the edge to the right
            MPI_Send(&nMigrate[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
            MPI_Send(&field.pos[0], nMigrate[1], MPI_DOUBLE, procID+1, sendMigrate, MPI_COMM_WORLD);

            // Receives the halo from the right
            MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
            MPI_Recv(&recvVectR, sizeRecvMigrate, MPI_DOUBLE, procID+1, recvMigrate, MPI_COMM_WORLD, status);
            // Removing previous particles
            resizeField(field, nMigrate[0]+nMigrate[1]);
            // Adding the particles to the field
            field.pos[0].insert(field.pos[0].end(), recvVectL.begin(), recvVectL.end());
            field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
        }
        else if(procID == nTasks - 1){ // No right neighbor domain
            // left halo: [l[0] , l[0] + boxSize[
            // left edge: [l[0] + boxSize , l[0] + 2*boxSize[
            // inner domain: [l[0] + 2*boxSize , u[0][
            computeMigrateIndex(field.pos[0],indexMigrate, nMigrate, l0+boxSize, u0);
            sortParticles(field, indexMigrate);

            if(procID%2==0){
                // Sends the edge to the left (procID is even)
                MPI_Send(&nMigrate[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[0], MPI_DOUBLE, procID-1, sendMigrate, MPI_COMM_WORLD);

                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvMigrate, MPI_DOUBLE, procID-1, recvMigrate, MPI_COMM_WORLD, status);
                // Removing previous particles
                resizeField(field, nMigrate[0]+nMigrate[1]);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].end(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }else{
                // Receives the edge from the left (procID is odd)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvMigrate, MPI_DOUBLE, procID-1, recvMigrate, MPI_COMM_WORLD, status);

                // Sends the edge to the left (procID is odd)
                MPI_Send(&nMigrate[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[0], MPI_DOUBLE, procID-1, sendMigrate, MPI_COMM_WORLD);
                // Removing previous particles
                resizeField(field, nMigrate[0]+nMigrate[1]);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].end(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }

        }
        else{ // Domain in the middle
            // left halo: [l[0] , l[0] + boxSize[
            // left edge: [l[0] + boxSize , l[0] + 2*boxSize[
            // inner domain: [l[0] + 2*boxSize , u[0] - 2*boxSize[
            // right edge: [u[0] - 2*boxSize , u[0] - boxSize[
            // right halo: [u[0] - boxSize , u[0]]
            computeMigrateIndex(field.pos[0],indexMigrate, nMigrate, l0+boxSize, u0-boxSize);
            sortParticles(field, indexMigrate);
            if(procID%2==0){
                // Sends the edge to the left (procID is even)
                MPI_Send(&nMigrate[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[0], MPI_DOUBLE, procID-1, sendMigrate, MPI_COMM_WORLD);
                // Sends the edge to the right (procID is even)
                MPI_Send(&nMigrate[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[1], MPI_DOUBLE, procID+1, sendMigrate, MPI_COMM_WORLD);

                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvMigrate, MPI_DOUBLE, procID-1, recvMigrate, MPI_COMM_WORLD, status);
                // Receives the edge from the right (procID is even)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectR, sizeRecvMigrate, MPI_DOUBLE, procID+1, recvMigrate, MPI_COMM_WORLD, status);
                // Removing previous particles
                resizeField(field, nMigrate[0]+nMigrate[1]);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].end(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }else{
                // Receives the edge from the right (procID is even)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID+1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectR, sizeRecvMigrate, MPI_DOUBLE, procID+1, recvMigrate, MPI_COMM_WORLD, status);
                // Receives the edge from the left (procID is even)
                MPI_Recv(&sizeRecvMigrate, 1, MPI_INT, procID-1, recvData, MPI_COMM_WORLD, status);
                MPI_Recv(&recvVectL, sizeRecvMigrate, MPI_DOUBLE, procID-1, recvMigrate, MPI_COMM_WORLD, status);
                // Sends the edge to the right (procID is even)
                MPI_Send(&nMigrate[1], 1, MPI_INT, procID+1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[1], MPI_DOUBLE, procID+1, sendMigrate, MPI_COMM_WORLD);
                // Sends the edge to the left (procID is even)
                MPI_Send(&nMigrate[0], 1, MPI_INT, procID-1, sendData, MPI_COMM_WORLD);
                MPI_Send(&field.pos[0], nMigrate[0], MPI_DOUBLE, procID-1, sendMigrate, MPI_COMM_WORLD);
                // Removing previous particles
                resizeField(field, nMigrate[0]+nMigrate[1]);
                // Adding the particles to the field
                field.pos[0].insert(field.pos[0].end(), recvVectL.begin(), recvVectL.end());
                field.pos[0].insert(field.pos[0].end(), recvVectR.begin(), recvVectR.end());
            }
        }
    }
    else{
        // Single processor
        return;
    }

}

void processUpdate(Field& localField, SubdomainInfo& subdomainInfo){
    // --- call sendMigrate ---
    /*
    shareMigrate(&localField, subdomainInfo.u[0], subdomainInfo.l[0],
         subdomainInfo.boxSize, subdomainInfo.procID, subdomainInfo.nTasks);
    // --- call shareOverlap ---
    shareOverlap(&localField, subdomainInfo.u[0], subdomainInfo.l[0],
         subdomainInfo.boxSize, subdomainInfo.procID, subdomainInfo.nTasks);
    */
}

void timeStepFinding(Field* currentField){
    // Find on each process the maximum acceptable time step

    // Reduce this information among all processes
}


void computeDomainIndex(std::vector<double> &posX,
    std::vector<double> &limits, std::vector<int> &nbPartNode,
    std::vector< std::pair<int,int> > &index, int nTasks){
    // Loop over particles
    for(int i=0 ; i<posX.size() ; i++){
        int domainNumber = getDomainNumber(posX[i], limits, nTasks);
        index.push_back(std::make_pair(domainNumber,i));
        nbPartNode[domainNumber]++;
    }
}

int getDomainNumber(double x, std::vector<double> &limits, int nTasks){
    // NAIVE !!! (implement a tree based search instead)
    int i;
    for(i=0 ; i<nTasks && x>limits[i] ; i++);
    return i-1;
}


// ------------------------------------------------------------
// FROM SORT.CPP

// left halo: [l[0] , l[0] + boxSize[
// left edge: [l[0] + boxSize , l[0] + 2*boxSize[
// inner domain: [l[0] + 2*boxSize , u[0] - 2*boxSize[
// right edge: [u[0] - 2*boxSize , u[0] - boxSize[
// right halo: [u[0] - boxSize , u[0]]
void computeMigrateIndex(std::vector<double>& posX,
    std::vector< std::pair<int,int> >& index, int* nMigrate,
    double Xmin, double Xmax){

    nMigrate[0]=0; nMigrate[1]=0;
    for(unsigned int i=0; i<posX.size(); ++i){
        if(posX[i]>Xmax){
            index.push_back( std::make_pair(rightMigrate,i) );
            ++nMigrate[1];
        }
        else if(posX[i]<Xmin){
            index.push_back( std::make_pair(leftMigrate,i) );
            ++nMigrate[0];
        }
        else{
            index.push_back( std::make_pair(noMigrate, i) );
        }
    }
}

void computeOverlapIndex(std::vector<double>& posX,
    std::vector< std::pair<int,int> >& index, int* nOverlap,
    double leftMinX, double leftMaxX, double rightMinX, double rightMaxX){

    for(unsigned int i=0; i<posX.size(); ++i){
        if(posX[i]>rightMinX && posX[i]<rightMaxX){ // Just in order to test, has to be adapted to the geometry
            index.push_back( std::make_pair(rightOverlap, i) );
            ++nOverlap[1];
        }
        else if(posX[i]<leftMaxX && posX[i]>leftMinX){
            index.push_back( std::make_pair(leftOverlap,i) );
            ++nOverlap[0];
        }
        else{
            index.push_back( std::make_pair(noOverlap, i) );
        }
    }
}

bool sortFunction(const std::pair<int,int>& one, const std::pair<int,int>& two){
    return (one.first < two.first);
}

void sortParticles(Field& field, std::vector< std::pair<int,int> >& index){
    // Sort the index vector
    std::sort(index.begin(), index.end(), sortFunction);
    int N=field.pos[0].size();
    // Temporary vectors for sorting
    std::vector<double> tmp(N);
    std::vector<int> tmpType(N);
    // --- Sorts all data one by one ---
    for(int coord=0 ; coord<3 ; coord++){
        // Position reordering
        for(int i=0; i<N; ++i)
            tmp[i]=field.pos[coord][ index[i].second ];
        (field.pos[coord]).swap(tmp);
        // Speed reordering
        for(int i=0; i<N; ++i)
            tmp[i]=field.speed[coord][ index[i].second ];
        (field.speed[coord]).swap(tmp);
    }
    // Density reordering
    for(int i=0; i<N; ++i)
        tmp[i]=field.density[ index[i].second ];
    (field.density).swap(tmp);
    // Pressure reordering
    for(int i=0; i<N; ++i)
        tmp[i]=field.pressure[ index[i].second ];
    (field.pressure).swap(tmp);
    // Mass reordering
    for(int i=0; i<N; ++i)
        tmp[i]=field.mass[ index[i].second ];
    (field.mass).swap(tmp);
    // Type reordering
    for(int i=0; i<N; ++i)
        tmpType[i]=field.type[ index[i].second ];
    (field.type).swap(tmpType);

    tmp.clear();// Ok?
    tmp.shrink_to_fit();// Ok?
    tmpType.clear();// Ok?
    tmpType.shrink_to_fit();// Ok?
}

void resizeField(Field& field, int nMigrate){
    int finalSize=(field.pos[0]).size()-nMigrate;
    for(int coord=0 ; coord<3 ; coord++){
        // Position resize
        (field.pos[coord]).resize(finalSize);
        (field.pos[coord]).shrink_to_fit();
        // Speed reordering
        (field.speed[coord]).resize(finalSize);
        (field.speed[coord]).shrink_to_fit();
    }
    // Density reordering
    (field.density).resize(finalSize);
    (field.density).shrink_to_fit();
    // Pressure reordering
    (field.pressure).resize(finalSize);
    (field.pressure).shrink_to_fit();
    // Mass reordering
    (field.mass).resize(finalSize);
    (field.mass).shrink_to_fit();
    // Type reordering
    (field.type).resize(finalSize);
    (field.type).shrink_to_fit();
}
