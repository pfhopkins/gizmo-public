/* file_io/hdf5_deflate_filter.cc — application-side registration of the HDF5
 * deflate (gzip) filter, backed by the linked zlib.
 *
 * Some HDF5 library builds are compiled WITHOUT zlib, so the builtin
 * H5Z_FILTER_DEFLATE filter is absent at runtime. Any deflate-compressed
 * dataset then fails to read — HDF5 reports "required filter 'deflate' is not
 * registered" — and (because H5Dread's return code was historically ignored)
 * GIZMO silently received zero/garbage data. GIZMO always links zlib (-lz),
 * so it can supply the filter itself.
 *
 * gizmo_register_hdf5_deflate_filter() installs a zlib-backed deflate filter
 * (filter id H5Z_FILTER_DEFLATE) ONLY when the HDF5 library does not already
 * provide one. On a normal HDF5 build (deflate present) it is a no-op, so
 * behavior is unchanged everywhere except on zlib-less HDF5 builds, where it
 * restores the ability to read/write compressed HDF5.
 *
 * Call once at startup (begrun), before any IC / snapshot read.
 *
 * Written for GIZMO.
 */
#include <hdf5.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "../declarations/allvars.h"
#include "../core/proto.h"

/* HDF5 filter callback: decompress on read (H5Z_FLAG_REVERSE), compress on
 * write. Contract: transform the `nbytes` bytes in *buf, free the old buffer,
 * store the new buffer + its allocation size, and return the valid byte
 * count (0 signals failure to HDF5). */
static size_t gizmo_zlib_deflate_filter(unsigned int flags, size_t cd_nelmts,
                                        const unsigned int cd_values[],
                                        size_t nbytes, size_t *buf_size, void **buf)
{
    if(flags & H5Z_FLAG_REVERSE)
    {
        /* ---- decompress (inflate) into a geometrically-grown buffer ---- */
        size_t nalloc = (nbytes ? nbytes : 1) * 4 + 16;
        void  *outbuf = malloc(nalloc);
        if(outbuf == NULL) return 0;

        z_stream strm;
        memset(&strm, 0, sizeof(strm));
        if(inflateInit(&strm) != Z_OK) { free(outbuf); return 0; }
        strm.next_in   = (Bytef *)(*buf);
        strm.avail_in  = (uInt)nbytes;
        strm.next_out  = (Bytef *)outbuf;
        strm.avail_out = (uInt)nalloc;

        for(;;)
        {
            int status = inflate(&strm, Z_SYNC_FLUSH);
            if(status == Z_STREAM_END) break;
            if(status != Z_OK) { inflateEnd(&strm); free(outbuf); return 0; }
            if(strm.avail_out == 0)
            {
                size_t used = nalloc;
                nalloc *= 2;
                void *nb = realloc(outbuf, nalloc);
                if(nb == NULL) { inflateEnd(&strm); free(outbuf); return 0; }
                outbuf = nb;
                strm.next_out  = (Bytef *)outbuf + used;
                strm.avail_out = (uInt)(nalloc - used);
            }
        }
        size_t outsize = (size_t)strm.total_out;
        inflateEnd(&strm);
        free(*buf);
        *buf = outbuf;
        *buf_size = nalloc;
        return outsize;
    }
    else
    {
        /* ---- compress (deflate) ---- */
        int level = (cd_nelmts > 0) ? (int)cd_values[0] : Z_DEFAULT_COMPRESSION;
        uLongf bound = compressBound((uLong)nbytes);
        void  *outbuf = malloc(bound ? (size_t)bound : 1);
        if(outbuf == NULL) return 0;
        if(compress2((Bytef *)outbuf, &bound, (const Bytef *)(*buf), (uLong)nbytes, level) != Z_OK)
            { free(outbuf); return 0; }
        free(*buf);
        *buf = outbuf;
        *buf_size = (size_t)bound;
        return (size_t)bound;
    }
}

void gizmo_register_hdf5_deflate_filter(void)
{
    /* Normal case: HDF5 was built with zlib and already provides deflate.
     * Do nothing — leave the builtin filter in place. */
    if(H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0) return;

    static const H5Z_class2_t gizmo_deflate_class = {
        H5Z_CLASS_T_VERS,                  /* H5Z_class_t struct version */
        (H5Z_filter_t)H5Z_FILTER_DEFLATE,  /* filter id (1 = deflate)    */
        1, 1,                              /* encoder + decoder present  */
        "gizmo-zlib-deflate",              /* filter name for diagnostics*/
        NULL, NULL,                        /* can_apply / set_local      */
        gizmo_zlib_deflate_filter          /* the filter operation       */
    };
    if(H5Zregister(&gizmo_deflate_class) < 0)
    {
        if(ThisTask == 0)
            printf("WARNING: HDF5 lacks the deflate filter and registering the "
                   "GIZMO zlib-backed replacement failed; compressed HDF5 files "
                   "will not be readable.\n");
        return;
    }
    if(ThisTask == 0)
        printf("HDF5 deflate filter absent from libhdf5; registered GIZMO "
               "zlib-backed deflate filter (compressed HDF5 I/O enabled).\n");
}
