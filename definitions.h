//
//  definitions.h
//
//  Created by Thomas Wong on 06/11/19.
//  Copyright © 2019 Thomas Wong. All rights reserved.
//

#ifndef definitions_h
#define definitions_h


#define COVER_THRES 5 // the minimum coverage for updating the reads

#define MAPQ_THRES 20 // the alignments with mapq < MAPQ_THRES are ignored

#define MAP_TRIMLEN_RATIO_THRES 0.5 // the alignment is ignored if more than MAP_TRIMLEN_RATIO_THRES of bases are trimed

#define QV_THRES 25 // min quality value of a base to consider

#define QUAL_VALUE_CHANGED_BASE 40 // the quality value of a base being updated

#endif
