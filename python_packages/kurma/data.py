import kurma.frames as fr

class trajectory:
    frames = None
    names = None
    box = None
    natoms = None
    nframes = None
    fn = None

    def __init__(self,inpf='',fn=[]):
        if inpf == '':
            pass
        else:
            if fn == []:
                self.nframes = fr.get_frame_total(inpf)
                fn = range(1,self.nframes+1)
            else:
                pass

            self.natoms = fr.get_natoms(inpf)
            self.frames,self.box,self.names = fr.get_frames(inpf,fn)
            self.fn = fn[0:len(self.frames)]
            self.nframes = len(self.frames)
