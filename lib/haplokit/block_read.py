"""
    block_reader module
    ~~~~~~~~~~~~~~~~~~~

    Identify LD blocks from haploview LD files.
"""

import os
import glob

from collections import defaultdict



class BlockIdentifier:
    """Indentify haplo-blocks automatically from haploview LD files.

    :param assoc_inst: an `Assoc` instance."""
    def __init__(self, assoc_inst):
        self.config = assoc_inst.config
        self.reportdir = os.path.join(self.config.get('ROUTINE', None), 'report')
        self.snpinfo = self.config.get('SNPFILE', None)
        self.dprime_cutoff = self.config.get('D_CUTOFF', None) or 0.9
        self.blocks = {}

    def block_go(self):
        LD_files = self.read_LDs()
        for f in LD_files:
            self.dprime_reader(f)
        return self.putdown()

    def read_LDs(self):
        """Iterate all `*LD` files in haploview dir."""
        haploview_dir = os.path.join(self.reportdir, 'haploview')
        return glob.glob('%s/*LD' % haploview_dir)

    def putdown(self):
        output = os.path.join(os.path.dirname(self.snpinfo), 'hap.txt')
        fmt = '**\t{0}\t{1}\n'
        with open(output, 'wt') as fh:
            for gene in self.blocks:
                blocks = self.blocks[gene]
                count = len(blocks)
                if count == 0:
                    continue
                elif count == 1:
                    fh.write(fmt.format(gene, '\t'.join(blocks[0])))
                else:
                    for n, b in enumerate(blocks):
                        name = gene + '-' + str(n + 1)
                        fh.write(fmt.format(name, '\t'.join(b)))
        return output

    def dprime_reader(self, filename):
        """Main function of this module."""
        genename = os.path.basename(filename).split('.')[0]
        raw_blocks = []
        dc = DprimeCollecter(filename, self.dprime_cutoff)
        dc.load_data()
        for site in dc.all_sites:
            #print('====site=====', site)
            try:
                # this is biggest possible block of this site.
                # we will narrow the range by checking if the 
                # sites in this block are linked.
                tmpblock = [site] + dc.dprimes[site]
                count = len(tmpblock)
                if count == 2:
                    # only 2 sites in the block, it must be true one.
                    raw_blocks.append(tmpblock)
                else:
                    # iterate the sites from left to right cause they're ordered.
                    # once it's breaked, the rest sites are given up and the former
                    # ones are count as a block.
                    n = 1
                    while n < count - 1:
                        next_one = tmpblock[n]
                        try:
                            lcs = [next_one] + dc.dprimes[next_one]
                            #if lcs == tmpblock[n:]:
                            #    n += 1
                            #    continue
                            if lcs == tmpblock[n:] or set(lcs) & set(tmpblock[n:]) == set(tmpblock[n:]):
                                # current site is linked with all the rest sites of th possible block
                                if n == count - 2:
                                    # reach last, all sites in the tmpblock are linked, it's a true block
                                    raw_blocks.append(tmpblock)
                                    break
                                else:
                                    # current site checked, move forward.
                                    n += 1
                                    continue
                            else:
                                # current site fail the check, tmpblock breaks here.
                                n += len(lcs)
                                raw_blocks.append(tmpblock[:n])
                                break
                        except KeyError:
                            # current site not linked with any ones, break here.
                            raw_blocks.append(tmpblock[:n + 1])
                            break
                        else:
                            raw_blocks.append(tmpblock)
            except KeyError:
                pass
        admit_blocks = self.block_merge(raw_blocks)
        self.blocks[genename] = admit_blocks

    @staticmethod
    def block_merge(blocks):
        """Since `dprime_reader` listed all possible blocks, it is very likely
        to have small blocks covered by bigger ones, the merging of these blocks
        is necessary.
        """
        count = len(blocks)
        if count == 1:
            return blocks
        else:
            merged = []
            point = 0
            mobile = 0
            while mobile < count or point < count:
                mobile = max(mobile, point + 1)
                ahead_block = blocks[point]
                try:
                    follow_block = blocks[mobile]
                    # blocks are odered by pisition, and followers is allways less
                    # then the former by at least one snv, so just need to consider
                    # this situation that follower covered by former.
                    if set(follow_block) & set(ahead_block) == set(follow_block):
                        mobile += 1
                        # here means that we reach the last block and it's contained
                        # in the ahead one, so not hesitate to push the ahead block 
                        # and break.
                        if mobile == count - 1:
                            merged.append(ahead_block)
                            break
                    else:
                        merged.append(ahead_block)
                        point = mobile
                except IndexError:
                    # point reaches last
                    if set(blocks[point]) & set(blocks[point - 1]) != set(blocks[point - 1]):
                        merged.append(blocks[point])
                    break
            return merged


class DprimeCollecter:
    """Extract information of position-ordered sites and d-prime values.

    :param LDfile: a LD file.
    :param cutoff: d-prime cutoff to decide blocks."""
    def __init__(self, LDfile, cutoff):
        self.filename = LDfile
        self.cutoff = cutoff
        self.dprimes = {}
        self.all_sites = []

    def load_data(self):
        count = 0
        record_all_sites = True
        record_ld_sites = True
        with open(self.filename, 'rt') as fh:
            for line in fh:
                count += 1
                if count == 1:
                    continue
                arr = line.split()
                if count == 2:
                    last_site = arr[0]
                    self.all_sites.append(arr[0])
                if arr[0] != last_site:
                    last_site = arr[0]
                    record_all_sites = False
                    record_ld_sites = True
                if record_all_sites:
                    self.all_sites.append(arr[1])

                # The largest block contains current sites reach here.
                # Because of the breakage, the rest sites won't form a
                # block with this sites, sou it's no need to record.
                if record_ld_sites:
                    if float(arr[2]) >= self.cutoff:
                        self.dprimes.setdefault(arr[0], []).append(arr[1])
                    else:
                        record_ld_sites = False
        return None



