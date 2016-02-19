import struct
import numpy
import re

#TODO: class

def load_Nanonics_NSOM_map(filename):
    channels = dict()
    with open(filename, 'rb') as f:
        data = f.read()
        type, header, data = re.split(b'-Start Header-|-End Header-', data)
        assert b"NAN File" in type

        # Parse Header
        header, comment = header.split(b',comment')  # "Comment" format is hard to parse, not touching it for now
        header = re.sub(b"[ \n\t\s\r]", b"", header)
        header = header.split(b',')
        header = [item.split(b'=') for item in header]
        header_dict = dict()
        header_dict[b'comment'] = comment
        for key, val in header:
            if b'.' in val:
                val = float(val)
            else:
                val = int(val)
            header_dict[key] = val

        #
        #Parse Channels
        #
        data = re.split(b'-Start Channel Header-', data)[1:]
        for dat in data:
            c_header, c_data = dat.split(b'-End Channel Header-')

            #c_header = re.sub(b"-Start Channel Header-", b"", c_header)
            c_header = re.sub(b"[ \n\t\s\r]", b"", c_header)

            c_header = c_header.split(b',')
            c_header = [item.split(b'=') for item in c_header]

            channel_header = dict()
            for key, val in c_header:
                if key == b"CMN" or key == b"CMX":
                    val = float(val)
                channel_header[key] = val

            channel_array = []

            fmt = '>H'
            size = struct.calcsize(fmt)
            bytelen = len(c_data)

            for j in range(int(bytelen / size)):
                byte = struct.unpack_from(fmt, c_data, size*j)[0]
                channel_array.append(byte)  # = struct.unpack_from("H", c_data)
            channel_size = len(channel_array)

            cmx = channel_header[b'CMX']
            cmn = channel_header[b'CMN']
            ref = header_dict[b'ReF']
            res = header_dict[b'ReS']
            a = (cmn-cmx)/(numpy.max(channel_array)-numpy.min(channel_array))
            b = cmx-a*numpy.max(channel_array)
            channel_array = numpy.multiply(channel_array, a) + b

            if channel_size < ref*2*res:
                channel_array = numpy.append(channel_array, cmn*numpy.ones(ref*2*res-channel_size))

            channel_array = numpy.reshape(channel_array, [res, ref*2])


            trace = channel_array[:, :ref]
            retrace = channel_array[:, ref+1:]

            channel = dict()
            channel['trace'] = trace
            channel['retrace'] = retrace
            channels[channel_header[b'CHN'].decode('ascii')] = channel

    return channels
