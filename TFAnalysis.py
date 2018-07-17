import numpy as np
import os
import subprocess
from tqdm import tqdm
from multiprocessing import Process, Manager, Pipe


# SeqPath = '/home/invites/jmorlot/HDD/Sequence/Sequence.npy'
# TFdatabasepath = '/home/invites/jmorlot/HDD/Datasets/TF_HOCOMOCOv10/TF_HOCOMOCO.txt'
# TFmotifpath = '/home/invites/jmorlot/HDD/Datasets/TF_HOCOMOCOv10/db/HUMAN/HOCOMOCOv10_HUMAN_mono_meme_format.meme'

class TF_enrichement():
    def __init__(self,SeqPath,TFdatabasepath,TFmotifpath):

        f = open(TFdatabasepath,'r')
        TF_database = [l.split('\n')[0] for l in f]
        f.close()
        self.TF_database = np.array(TF_database)

        self.Sequence = np.load(SeqPath)

        self.TFmotifpath = TFmotifpath

        if not os.path.exists('TF_Analysis/'): os.mkdir('TF_Analysis/')

    def ComputeTFEnrichement(self,matrix,idxNull,NR=30):
        '''
            Compute the enrichement of TF compared to a null model
            matrix: Binary matrix where each lines correspond to the position to check
            idxNull: Position of all possible peaks
            NR: Number of null models
        '''

        print 'Get TF'
        indexes = [np.where(matrix[i,:]>0)[0] for i in range(matrix.shape[0])]
        # self.NI = len(indexes)

        self.TF = self.GetFimoTF(indexes)

        print 'Null model'
        self.TFR = np.zeros((len(self.TF_database),NR,len(indexes)))
        # self.TFR = np.zeros((len(indexes),len(self.TF_database),NR))
        # for nr in tqdm(range(NR)):
        for i,idx in tqdm(enumerate(indexes)):
            #Generate false indexes conserving the number of peaks
            indexesR = [idxNull[np.random.permutation(len(idxNull))[:len(idx)]] for nr in range(NR)]
            self.TFR[:,:,i] = self.GetFimoTF(indexesR)

        #Zscore
        m = self.TFR.mean(axis=1)
        s = self.TFR.std(axis=1)
        s[s==0] = 1
        s[s==0] = 1
        self.TFZ = (self.TF - m)/s

        return self.TFZ

    def GetFimoTF(self,indexes):
        print '\tWrite sequences in a FASTA File'
        self.seqpaths = self.WriteSequenceFile(indexes)

        print '\tLaunch FIMO in parallel on sequences'
        self.FIMO_MP(self.seqpaths)

        print '\tRead FIMO file and construct TF matrix'
        self.fimopaths = [seqpath+'fimo.txt' for seqpath in self.seqpaths]
        TF = self.ReadFimoFile(self.fimopaths)

        return TF

    def WriteSequenceFile(self,indexes):
        # self.seqnum = []
        seqpaths = []
        #print 'Write Sequence File'
        for k,idx in tqdm(enumerate(indexes)):
            seq = self.Sequence[idx]

            seqpath = 'TF_Analysis/sequence{}/'.format(k)
            if not os.path.exists(seqpath): os.mkdir(seqpath)
            f = open(seqpath+'sequences.txt','w')
            for i,n in enumerate(seq):
                f.write('>seq' + str(i) + '\n')
                f.write(n + '\n\n')
            f.close()
            seqpaths.append(seqpath)

        return seqpaths

    def FIMO_MP(self,seqpaths):
        processes = []
        for seqpath in seqpaths:
            p = Process(target=launchFIMO, args=(self.TFmotifpath,seqpath))
            processes.append(p)
            p.start()

        for process in processes:
            process.join()

    def ReadFimoFile(self,fimopaths):
        '''
            Get a list of TF from a list of genomic positions usinf FIMO
        '''

        TF = np.zeros((len(self.TF_database),len(fimopaths)))

        #print 'Get FIMO result'
        for k,fimopath in tqdm(enumerate(fimopaths)):
            TF_idx = []
            print k
            # seqnum_i = []
            f = open(fimopath,'r')
            for l in f:
                try:
                    tf = ((l.split('\n')[0]).split('\t')[0]).split('_HUMAN')[0]
                    tf = np.where(self.TF_database==tf)[0][0]
                    # seqn = int((l.split('\t')[2]).split('seq')[1])
                    TF_idx.append(tf)
                    # seqnum_i.append(seqn)
                except:
                    pass

            f.close()

            TF_idx = np.array(TF_idx)
            TFU,counts = np.unique(TF_idx,return_counts=True)
            TF[TFU,k] = counts

        return TF

def launchFIMO(TFmotifpath,SequencesPath):
    '''
        Launch FIMO using the motif database, loacated at TFmotifpath,
        on the sequences located at SequencesPath in FASTA format
    '''
    CMD = 'fimo -oc '+SequencesPath+' ' + TFmotifpath + ' ' + SequencesPath+'sequences.txt'
    subprocess.call(CMD,shell=True,stderr=subprocess.STDOUT)
