
from pt_amino_acids import amino_acids
from pt_tcr_distances_blosum import blosum, bsd4

class DistanceParams:
    def __init__(self, config_string=None ):
        self.gap_penalty_v_region = 4
        self.gap_penalty_cdr3_region = 8
        self.weight_v_region = 1
        self.weight_cdr3_region = 3
        self.distance_matrix = bsd4
        self.align_cdr3s = False
        self.trim_cdr3s = True
        self.scale_factor = 1.0
        if config_string:
            l = config_string.split(',')
            for tag,val in [x.split(':') for x in l ]:
                if tag == 'gap_penalty_cdr3_region':
                    self.gap_penalty_cdr3_region = float(val)
                elif tag == 'gap_penalty_v_region':
                    self.gap_penalty_v_region = float(val)
                elif tag == 'weight_cdr3_region':
                    self.weight_cdr3_region = float(val)
                elif tag == 'weight_v_region':
                    self.weight_v_region = float(val)
                elif tag == 'scale_factor':
                    self.scale_factor = float(val)
                elif tag == 'align_cdr3s':
                    assert val in ['True','False']
                    self.align_cdr3s = ( val == 'True' )
                elif tag == 'trim_cdr3s':
                    assert val in ['True','False']
                    self.trim_cdr3s = ( val == 'True' )
                else:
                    print('unrecognized tag:',tag)
                    assert False
            print('config_string: {} self: {}'.format( config_string, self ))

    def __str__(self):
        return 'DistanceParams: gap_penalty_v_region= {} gap_penalty_cdr3_region= {} weight_v_region= {} weight_cdr3_region= {} align_cdr3s= {} trim_cdr3s= {}'\
            .format( self.gap_penalty_v_region, self.gap_penalty_cdr3_region,
                     self.weight_v_region, self.weight_cdr3_region,
                     self.align_cdr3s, self.trim_cdr3s )

default_distance_params = DistanceParams()

gap_character = "."

def blosum_character_distance( a, b, gap_penalty, params ):
    if a== gap_character and b == gap_character:
        return 0
    elif a == '*' and b == '*':
        return 0
    elif a == gap_character or b == gap_character or a=='*' or b=='*':
        return gap_penalty
    else:
        # assert a in amino_acids
        # assert b in amino_acids
        # maxval = min( blosum[(a,a)], blosum[(b,b)] )
        # return maxval - blosum[(a,b)]
        return params.distance_matrix[ (a,b) ]

def blosum_sequence_distance( aseq, bseq, gap_penalty, params ):
    assert len(aseq) == len(bseq)
    dist = 0.0
    for a,b in zip(aseq,bseq):
        if a == ' ':
            assert b== ' '
        else:
            dist += blosum_character_distance( a, b, gap_penalty, params )
    return dist

def align_cdr3s( a, b, gap_character ):
    if len(a) == len(b):
        return (a[:],b[:])

    if len(a)<len(b): ## s0 is the shorter sequence
        s0,s1 = a,b
    else:
        s0,s1 = b,a

    lendiff = len(s1)-len(s0)

    best_score=-1000

    # the gap comes after s0[gappos]
    for gappos in range(len(s0)-1):
        score=0
        for i in range(gappos+1):
            score += blosum[ (s0[i],s1[i]) ]
        for i in range(gappos+1,len(s0)):
            score += blosum[ (s0[i],s1[i+lendiff]) ]
        if score>best_score:
            best_score = score
            best_gappos = gappos
    ## insert the gap
    s0 = s0[:best_gappos+1] + gap_character*lendiff + s0[best_gappos+1:]

    assert len(s0) == len(s1)

    if len(a)<len(b): ## s0 is the shorter sequence
        return ( s0, s1 )
    else:
        return ( s1, s0 )

## align
##
##   shortseq[        ntrim: gappos   ] with longseq[ ntrim: gappos ] and
##   shortseq[ -1*remainder: -1*ctrim ] with longseq[ -1*remainder: -1*ctrim ]
##
## but be careful about negative indexing if ctrim is 0
##
## the gap comes after position (gappos-1) ie there are gappos amino acids before the gap
##
##
## DOES NOT INCLUDE THE GAP PENALTY
##
def sequence_distance_with_gappos( shortseq, longseq, gappos, params ):
    ntrim = 3 if params.trim_cdr3s else 0
    ctrim = 2 if params.trim_cdr3s else 0
    remainder = len(shortseq)-gappos
    dist = 0.0
    count =0
    if ntrim < gappos:
        for i in range(ntrim,gappos):
            dist += params.distance_matrix[ (shortseq[i], longseq[i] ) ]
            count += 1
    if ctrim < remainder:
        for i in range(ctrim, remainder):
            dist += params.distance_matrix[ (shortseq[-1-i], longseq[-1-i] ) ]
            count += 1
    return dist,count


def weighted_cdr3_distance( seq1, seq2, params ):
    shortseq,longseq = (seq1,seq2) if len(seq1)<=len(seq2) else (seq2,seq1)

    ## try different positions of the gap
    lenshort = len(shortseq)
    lenlong = len(longseq)
    lendiff = lenlong - lenshort

#    assert lenshort>3 ##JCC testing
    assert lenshort > 1##JCC testing
    assert lendiff>=0
    if params.trim_cdr3s:
        assert lenshort > 3+2 ## something to align...

    if not params.align_cdr3s:
        ## try to replicate old (strange) behavior: "gap_spot = min( 3, len(shortseq)/2 )" in ../tcr_distances.py
        ## shortseq in that code had already been trimmed by 3,2 residues
        gappos = min( 6, 3 + (lenshort-5)//2 )
        best_dist,count = sequence_distance_with_gappos( shortseq, longseq, gappos, params )

    else:
        ## the CYS and the first G of the GXG are 'aligned' in the beta sheet
        ## the alignment seems to continue through roughly CYS+4
        ## ie it's hard to see how we could have an 'insertion' within that region
        ## gappos=1 would be a insertion after CYS
        ## gappos=5 would be a insertion after CYS+4 (5 rsds before the gap)
        ## the full cdr3 ends at the position before the first G
        ## so gappos of len(shortseq)-1 would be gap right before the 'G'
        ## shifting this back by 4 would be analogous to what we do on the other strand, ie len(shortseq)-1-4
        min_gappos = 5
        max_gappos = len(shortseq)-1-4
        while min_gappos>max_gappos:
            min_gappos -= 1
            max_gappos += 1
        for gappos in range( min_gappos, max_gappos+1 ):
            dist, count = sequence_distance_with_gappos( shortseq, longseq, gappos, params )
            if gappos>min_gappos:
                assert count==best_count
            if gappos == min_gappos or dist < best_dist:
                best_dist = dist
                best_gappos = gappos
                best_count = count

    return  params.weight_cdr3_region * best_dist + lendiff * params.gap_penalty_cdr3_region


def compute_distance(t1,t2,chains,rep_dists,distance_params): #    t1/2 = [ va_reps, vb_reps, l['cdr3a'], l['cdr3b'] ]
    dist=0.0
    if 'A' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[0] for y in t2[0] ) ) +\
                weighted_cdr3_distance( t1[2], t2[2], distance_params )
    if 'B' in chains:
        dist += min( ( rep_dists[x][y] for x in t1[1] for y in t2[1] ) ) +\
                weighted_cdr3_distance( t1[3], t2[3], distance_params )
    return distance_params.scale_factor * dist

