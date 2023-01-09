import numpy as np
from ..io import is_notebook

from tqdm import *
#from tqdm.notebook import tqdm as ipython_tqdm

#from tqdm.autonotebook import tqdm

class SASA_Results():
    def __init__(self, is_a_frame):
        self.resname  = None
        self.resID    = None
        self.chainID  = None
        self.sasa     = None
        self.nframe   = 0
        self.nresidue = 0
        self.raw = None
        self.is_a_frame = is_a_frame
    # Add a progress bar here too, it takes time
    def parse(self, data):
        if self.is_a_frame:
            self._parse_a_frame(data)
        else : 
            self._parse_raw(data)
    def _parse_a_frame(self, sasa_multi_frame_data):
         #bar_constructor = tqdm if not is_notebook() else ipython_tqdm
        # Above line is not compatible with ipywidet >= 8.0 for now ...
        bar_constructor = tqdm
        c = 0
        for _ in sasa_multi_frame_data:
            c += len(_['sasa'])
        print(f"Compiling results over {c} snapshots")
        with bar_constructor(total=self.nframe) as pbar:
            for result in sasa_multi_frame_data:
                for k in ['resname', 'resID', 'chainID']:
                    if getattr(self, k) is None:
                        setattr(self, k, np.array(result[k]))
                for curr_pose_sasa_list in result['sasa']:
                    self.nframe += 1
                    self.nresidue = len(curr_pose_sasa_list)
                    self.sasa = np.array(curr_pose_sasa_list) if self.sasa is None\
                                        else np.vstack( [self.sasa, curr_pose_sasa_list] )
                    pbar.update(1)
    def _parse_raw(self, list_sasa_dict):
        self.raw = list_sasa_dict
    ## Add a progres bar if necessary but not sure
    def stats(self):
        if not self.is_a_frame:
            raise TypeError("Cannot process statistics on non frame results, still ou can inspect me through \".raw\" attribute")

        stats = [ ( i_tup, np.zeros(self.nframe) )\
            for i_res, i_tup in enumerate(zip(self.resname, self.resID, self.chainID) )\
        ]

        for iframe in range(0, self.nframe):
            for ires in range(0, self.nresidue):
                cur_index = iframe * self.nresidue + ires
                (cur_asa_abs, cur_asa_frac) = self.sasa[cur_index]
                stats[ires][1][iframe] = cur_asa_abs
                
        values = [ (i_res_tup, np.mean(data), np.std(data)) \
            for (i_res_tup, data) in stats ]
        
        return values
            