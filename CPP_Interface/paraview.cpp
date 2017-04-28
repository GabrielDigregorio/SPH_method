#include "paraview.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <stdint.h>
#include "paraview.h"
#include "swapbytes.h"

#ifdef USE_ZLIB
#include <zlib.h>
#else
#define Z_OK 0
#endif

const int __one__ = 1;
const bool isCpuLittleEndian = 1 == *(char*)(&__one__); // CPU endianness

// this routine writes a vector of double as a vector of float to a legacy VTK file f.
//  the vector is converted to float (32bits) and to "big endian" format (required by the legacy VTK format)

void write_vectorLEGACY(std::ofstream &f, 
                        std::vector<double> const * pos, int dim, int nbpStart, int nbpEnd, bool binary)
{
    /*std::cout << "write_vectorLEGACY";
    std::cout << "dim=" << dim << '\n';
    std::cout << "nbpStart=" << nbpStart << '\n';
    std::cout << "nbpEnd=" << nbpEnd << '\n';*/
    if(!binary) 
    {
        for(int i=nbpStart; i<nbpEnd; ++i)
        {
            for(int j=0;j<dim;++j)
                f << pos[j][i] << " ";
            f << '\n';
        }
    }
    else
    {
        if(isCpuLittleEndian)
            for(int i=nbpStart; i<nbpEnd; ++i)
            {
                // float+little endian => double should be converted to float, then swapped
                for(int j=0; j<dim; ++j)
                {
                    float fx = (float)pos[j][i]; 
                    uint32_t x = swap_uint32(*(uint32_t*)&fx); // convert if CPU is little endian
                    f.write((char*)&x, sizeof(uint32_t));
                }
            }
        else
        {
            // double+bigendian => vector can be written as in memory
            //f.write(reinterpret_cast<char const*>(&(pos[0])), pos.size()*sizeof(double));
            // float+bigendian => vector should be converted to float
            for(int i=nbpStart; i<nbpEnd; ++i)
                for(int j=0; j<3; ++j)
                {
                    float fx = (float)pos[j][i]; 
                    f.write((char*)&fx, sizeof(uint32_t));
                }
        }
    }
    //std::cout << "done.\n";
}

// export results to paraview (VTK polydata - legacy file fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]
//   binary:   'true' for binary format, 'false' for ASCII

void paraviewLEGACY(std::string const &filename, 
              int step,
              std::vector<double> const (&pos)[3],
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> (*)[3]> const &vectors,
              int nbpStart, int nbpEnd,  
              bool binary)
{
    int nbp = nbpEnd-nbpStart;
    
    // build file name + stepno + vtk extension
    std::stringstream s; s << "Results/" << filename << "_" << std::setw(8) << std::setfill('0') << step << ".vtk";

    // open file
    //std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    f << std::scientific;
    // header
    f << "# vtk DataFile Version 3.0\n";
    f << "file written by sph.exe\n";
    f << (binary ? "BINARY\n" : "ASCII\n");
    f << "DATASET POLYDATA\n";

    // points
    f << "POINTS " << nbp << " float\n";
    write_vectorLEGACY(f, pos, 3, nbpStart, nbpEnd, binary);
    
    // vertices
    f << "VERTICES " << nbp << " " << 2*nbp << "\n";
    if(!binary) 
    {    
        for(int i=0; i<nbp; ++i)
            f << "1 " << i-nbpStart << '\n';
        f << '\n'; // empty line (required)
    }
    else
    {
        int32_t type = isCpuLittleEndian? swap_int32(1) : 1;
        for(int i=0; i<nbp; ++i)
        {
            int32_t ii = isCpuLittleEndian? swap_int32(i) : i;
            f.write((char*)&type, sizeof(int));
            f.write((char*)&ii, sizeof(int));
        }
    }

    // fields
    f << "POINT_DATA " << nbp << '\n';
    f << "FIELD FieldData " << scalars.size()+vectors.size() << '\n';

    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        //assert(it->second->size()==nbp);
        f << it->first << " 1 " << nbp << " float\n";
        write_vectorLEGACY(f, &*it->second, 1, nbpStart, nbpEnd, binary);
    }

    // vector fields
    std::map<std::string, std::vector<double> (*)[3]>::const_iterator itV = vectors.begin();
    for(; itV!=vectors.end(); ++itV)
    {
        //assert(it->second->size()==3*nbp);
        f << itV->first << " 3 " << nbp << " float\n";
        write_vectorLEGACY(f, &(*itV->second)[0], 3, nbpStart, nbpEnd, binary);
    }
    f.close();
}


std::string zlibstatus(int status)
{
#ifdef USE_ZLIB
    switch(status)
    {
        case Z_OK: return "Z_OK";
        case Z_BUF_ERROR: return "Z_BUF_ERROR";
        case Z_MEM_ERROR: return "Z_MEM_ERROR";
        case Z_STREAM_ERROR: return "Z_STREAM_ERROR";
        default:
            std::stringstream str; str << "Unknown ("<< status << ")";
            return str.str();
    }
#else
    return "zlib missing";
#endif
}


size_t write_vectorXML(std::ofstream &f, std::vector<double> const * pos, int dim,
                       int nbpStart, int nbpEnd, bool usez)
{
    /*std::cout << "write_vectorXML\n";
    std::cout << "dim="<<dim << '\n';
    std::cout << "nbpStart="<<nbpStart << '\n';
    std::cout << "nbpEnd="<<nbpEnd << '\n';*/

    size_t written=0;
    int nbp = nbpEnd-nbpStart;

    if(!usez)
    {
        // data block size
        uint32_t sz = nbp*dim*sizeof(float);
        f.write((char*)&sz, sizeof(uint32_t)); written+=sizeof(uint32_t);
        // data
        for(int i=nbpStart; i<nbpEnd; ++i)
        {
            for(int j=0;j<dim; ++j)
            {
                float fx = (float)pos[j][i]; 
                f.write((char*)&fx, sizeof(float));
            }
        }
        written+=sz;
    }
    else
    {
        // convert double to float
        std::vector<float> buffer(nbp*dim);
        for(int i=nbpStart; i<nbpEnd; ++i)
            for(int j=0;j<dim; ++j)
                buffer[(i-nbpStart)*dim+j] = (float)pos[j][i];

        size_t sourcelen = buffer.size() *sizeof(float);
        size_t destlen = size_t(sourcelen * 1.001) + 12;  // see doc
        char *destbuffer = new char[destlen];
#ifdef USE_ZLIB
        int status = compress2((Bytef*)destbuffer, &destlen, 
                               (Bytef *)&(buffer[0]), sourcelen, Z_DEFAULT_COMPRESSION);
#else
        int status = Z_OK+1;
#endif
        if(status!=Z_OK)
        {
            std::cout << "ERROR: zlib Error status=" << zlibstatus(status) << "\n";
        }
        else
        {
            //std::cout << "block of size " << sourcelen << " compressed to " << destlen << '\n';
            // blocks description
            uint32_t nblocks=1;
            f.write((char*)&nblocks, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t srclen = (uint32_t)sourcelen;
            f.write((char*)&srclen, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t lastblocklen = 0;
            f.write((char*)&lastblocklen, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t szblocki = (uint32_t)destlen;
            f.write((char*)&szblocki, sizeof(uint32_t)); written += sizeof(uint32_t);
            // data
            f.write(destbuffer, destlen);
            written += destlen;
        }

        delete []destbuffer;
    }

    return written;
}


size_t write_vectorXML(std::ofstream &f, std::vector<int> const &pos, bool usez)
{
    //std::cout << "write_vectorXML(int)\n";
    size_t written=0;

    if(!usez)
    {
        // data block size
        uint32_t sz = pos.size()*sizeof(int);
        f.write((char*)&sz, sizeof(uint32_t)); written+=sizeof(uint32_t);
        // data
        for(int i=0; i<pos.size(); ++i)
        {
            int fx = pos[i]; 
            f.write((char*)&fx, sizeof(int));
        }
        written+=sz;
    }
    else
    {
        size_t sourcelen = pos.size() *sizeof(int);
        size_t destlen = size_t(sourcelen * 1.001) + 12;  // see doc
        char *destbuffer = new char[destlen];
#ifdef USE_ZLIB
        int status = compress2((Bytef *)destbuffer, &destlen, 
                               (Bytef *)&(pos[0]), sourcelen, Z_DEFAULT_COMPRESSION);
#else
        int status = Z_OK+1;
#endif
        if(status!=Z_OK)
        {
            std::cout << "ERROR: zlib Error status=" << zlibstatus(status) << "\n";
        }
        else
        {
            // blocks description
            uint32_t nblocks=1;
            f.write((char*)&nblocks, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t srclen = (uint32_t)sourcelen;
            f.write((char*)&srclen, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t lastblocklen = 0;
            f.write((char*)&lastblocklen, sizeof(uint32_t)); written += sizeof(uint32_t);
            uint32_t szblocki = (uint32_t)destlen;
            f.write((char*)&szblocki, sizeof(uint32_t)); written += sizeof(uint32_t);
            // data
            f.write(destbuffer, destlen);
            written += destlen;
        }

        delete []destbuffer;
    }

    return written;
}
 

// export results to paraview (VTK polydata - XML fomat)
//   filename: file name without vtk extension
//   pos:     positions (vector of size 3*number of particles)
//   step:    time step number
//   scalars: scalar fields defined on particles (map linking [field name] <=> [vector of results v1, v2, v3, v4, ...]
//   vectors: vector fields defined on particles (map linking [field name] <=> [vector of results v1x, v1y, v1z, v2x, v2y, ...]

// see http://www.vtk.org/Wiki/VTK_XML_Formats

void paraviewXML(std::string const &filename, 
                 int step,
                 std::vector<double> const (&pos)[3],
                 std::map<std::string, std::vector<double> *> const &scalars,
                 std::map<std::string, std::vector<double> (*)[3]> const &vectors,
                 int nbpStart, int nbpEnd,  
                 bool binary, 
                 bool usez)
{
#if !defined(USE_ZLIB)
    if(binary && usez)
    {
        //std::cout << "INFO: zlib not present - vtk file will not be compressed!\n";
        usez=false; 
    }
#endif

    int nbp = nbpEnd-nbpStart;
    
    // build file name + stepno + vtk extension
    std::stringstream s; s << "Results/" << filename << "_" << std::setw(8) << std::setfill('0') << step << ".vtp";
    std::stringstream s2; s2 << "Results/" << filename << "_" << std::setw(8) << std::setfill('0') << step << ".vtp.tmp";

    // open file
    //std::cout << "writing results to " << s.str() << '\n';
    std::ofstream f(s.str().c_str(), std::ios::binary | std::ios::out);
    std::ofstream f2(s2.str().c_str(), std::ios::binary | std::ios::out); // temp binary file
    f << std::scientific;

    size_t offset = 0;
    // header
    f << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"";
    f << ( isCpuLittleEndian? "LittleEndian" : "BigEndian") << "\" ";
    f << "header_type=\"UInt32\" "; // UInt64 should be better
    if(usez)
        f << "compressor=\"vtkZLibDataCompressor\" ";
    f << ">\n";
    f << "  <PolyData>\n";
    f << "    <Piece NumberOfPoints=\"" << nbp << "\" ";
    f << "NumberOfVerts=\"" << nbp << "\" ";
    f << "NumberOfLines=\"0\" ";
    f << "NumberOfStrips=\"0\" ";
    f << "NumberOfPolys=\"0\">\n";

    // ------------------------------------------------------------------------------------
    f << "      <PointData>\n";
    // scalar fields
    std::map<std::string, std::vector<double> *>::const_iterator it=scalars.begin();
    for(; it!=scalars.end(); ++it)
    {
        //assert(it->second->size()==nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << it->first << "\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, &*it->second, 1, nbpStart, nbpEnd, usez);
    }
    // vector fields
    std::map<std::string, std::vector<double> (*)[3]>::const_iterator itV = vectors.begin();
    for(; itV!=vectors.end(); ++itV)
    {
        //assert(it->second->size()==3*nbp);
        f << "        <DataArray type=\"Float32\" ";
        f << " Name=\"" << itV->first << "\" ";
        f << " NumberOfComponents=\"3\" ";
        f << " format=\"appended\" ";
        f << " RangeMin=\"0\" ";
        f << " RangeMax=\"1\" ";
        f << " offset=\"" << offset << "\" />\n";
        offset += write_vectorXML(f2, &(*itV->second)[0], 3, nbpStart, nbpEnd, usez);
    }
    f << "      </PointData>\n";

    // ------------------------------------------------------------------------------------
    f << "      <CellData>\n";
    f << "      </CellData>\n";  

    // ------------------------------------------------------------------------------------
    f << "      <Points>\n";
    f << "        <DataArray type=\"Float32\" ";
    f << " Name=\"Points\" ";
    f << " NumberOfComponents=\"3\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";   
    offset += write_vectorXML(f2, pos, 3, nbpStart, nbpEnd, usez);
    f << "      </Points>\n";
    // ------------------------------------------------------------------------------------
    f << "      <Verts>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"" << nbp-1 << "\" ";
    f << " offset=\"" << offset << "\" />\n";

    std::vector<int> connectivity(nbp);   // <= hard to avoid if zlib is used
    for(int i=0; i<nbp; ++i) connectivity[i]=i;
    offset += write_vectorXML(f2, connectivity, usez);

    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"1\" ";
    f << " RangeMax=\"" << nbp << "\" ";
    f << " offset=\"" << offset << "\" />\n";
    
    // reuse "connectivity" for offsets
    for(int i=0; i<nbp; ++i) connectivity[i]=i+1;
    offset += write_vectorXML(f2, connectivity, usez);

    f << "      </Verts>\n";

    std::vector<double> empty;
    // ------------------------------------------------------------------------------------
    f << "      <Lines>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez);
      
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez);
    f << "      </Lines>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Strips>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";   
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez);
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n"; 
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez);
    f << "      </Strips>\n";

    // ------------------------------------------------------------------------------------
    f << "      <Polys>\n";
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"connectivity\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";  
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez); 
    f << "        <DataArray type=\"Int32\" ";
    f << " Name=\"offsets\" ";
    f << " format=\"appended\" ";
    f << " RangeMin=\"0\" ";
    f << " RangeMax=\"1\" ";
    f << " offset=\"" << offset << "\" />\n";
    offset += write_vectorXML(f2, &empty, 1, 0, 0, usez); 
    f << "      </Polys>\n";

    f2.close();

    // ------------------------------------------------------------------------------------
    f << "    </Piece>\n";
    f << "  </PolyData>\n";
    // ------------------------------------------------------------------------------------
    f << "  <AppendedData encoding=\"raw\">\n";
    f << "    _";

    // copy temp binary file
    std::ifstream f3(s2.str().c_str(), std::ios::binary | std::ios::in);
    f << f3.rdbuf();
    f3.close();
    // remove temp file
    std::remove(s2.str().c_str());

    f << "  </AppendedData>\n";
    f << "</VTKFile>\n";

    f.close();
}

// interface

void paraview(std::string const &filename, 
              int step,
              std::vector<double> const (&pos)[3],
              std::map<std::string, std::vector<double> *> const &scalars,
              std::map<std::string, std::vector<double> (*)[3]> const &vectors,
              int nbpStart, int nbpEnd, 
              PFormat format)       
{
    switch(format)
    {
        case LEGACY_TXT:
            paraviewLEGACY(filename, step, pos, scalars, vectors, nbpStart, nbpEnd, false); break;
        case XML_BIN:
            paraviewXML(filename, step, pos, scalars, vectors, nbpStart, nbpEnd, true, false); break;
        case XML_BINZ:
            paraviewXML(filename, step, pos, scalars, vectors, nbpStart, nbpEnd, true, true); break;
        case LEGACY_BIN:
        default:
            paraviewLEGACY(filename, step, pos, scalars, vectors, nbpStart, nbpEnd, true); break;
    }
}      
