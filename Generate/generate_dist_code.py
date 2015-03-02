#!/usr/bin/env python
import StringIO,sys


def print_commons(line,IOobj):
    if len(line) < 72 or line[0] in ('c','C'):
        print >> IOobj, line
    else:
        split = line.split(',')
        end = split.pop()
        print >> IOobj, ','.join(split)
        print >> IOobj, '     & {}'.format(end)
def print_code(line,IOobj):
    if len(line) < 72 or line[0] in ('c','C'):
        print >> IOobj, line
    else:
        split = line[:72]
        end = line[72:]
        print >> IOobj, split
        print >> IOobj, '     & {}'.format(end)

distnames = ['cost3','cost4','cost3col','cost4col',
             'Y3','Y4','Y3col','Y4col',
             'eta3','eta4','eta3col','eta4col',
             'Yff','Pzff','rmass','beta','delY','Et','pT']
description = {
    'cost3':'anti-fermion CM frame polar angle',
    'cost4':'fermion CM frame polar angle',
    'cost3col':'anti-fermion lab frame polar angle',
    'cost4col':'fermion lab frame polar angle',
    'Y3':'anti-fermion CM frame rapidity',
    'Y4':'fermion CM frame rapidity',
    'Y3col':'anti-fermion lab frame rapidity',
    'Y4col':'fermion lab frame rapidity',
    'eta3':'anti-fermion CM frame pseudorapidity',
    'eta4':'fermion CM frame pseudorapidity',
    'eta3col':'anti-fermion lab frame pseudorapidity',
    'eta4col':'fermion lab frame pseudorapidity',
    'Yff':'ffbar system rapidity',
    'Pzff':'ffbar system z-momentum',
    'rmass':'ffbar invariant mass',
    'beta':'Lorentz boost',
    'delY':'ffbar rapidity difference',
    'Et':'fermion energy',
    'pT':'fermion tranverse momentum'
}
words = {
    'cost3':('cos3cm','\',f,\'bar (CM) Polar Angle','cos{/Symbol q}'),
    'cost4':('cos4cm','\',f,\' (CM) Polar Angle','cos{/Symbol q}'),
    'cost3col':('cos3','\',f,\'bar Polar Angle','cos{/Symbol q}'),
    'cost4col':('cos4','\',f,\' Polar Angle','cos{/Symbol q}'),
    'Y3':('y3cm','\',f,\'bar (CM) Rapidity','y'),
    'Y4':('y4cm','\',f,\' (CM) Rapidity','y'),
    'Y3col':('y3','\',f,\'bar Rapidity','y'),
    'Y4col':('y4','\',f,\' Rapidity','y'),
    'eta3':('eta3cm','\',f,\'bar (CM) Pseudorapidity','{/Symbol h}'),
    'eta4':('eta4cm','\',f,\' (CM) Pseudorapidity','{/Symbol h}'),
    'eta3col':('eta3','\',f,\'bar Pseudorapidity','{/Symbol h}'),
    'eta4col':('eta4','\',f,\' Pseudorapidity','{/Symbol h}'),
    'Yff':('y\',ff,\'','y_{\',ff,\'}','y_{\',ff,\'}'),
    'Pzff':('pz\',ff,\'','p^{z}_{\',ff,\'}','p^{z}_{\',ff,\'}'),
    'rmass':('M\',ff,\'','\',mff,\'','\',mff,\''),
    'beta':('beta','{\Symbol b}','{\Symbol b}'),
    'delY':('dely','{\Symbol D}y','{\Symbol D}y'),
    'Et':('E\',f,\'','E_\',f,\'','E_\',f,\''),
    'pT':('pT','p_T^\',f,\'','p_T')
}




inc_files = []
setting = StringIO.StringIO()
for name in distnames:
    common, binning, plotting, filling  = StringIO.StringIO(), StringIO.StringIO(), StringIO.StringIO(), StringIO.StringIO()
    
    print_commons('C {}'.format(description[name]),common)
    print_commons('      real*8 {0}max,{0}min,{0}w'.format(name),common)
    print_commons('      common/ext_{0}/{0}max,{0}min,{0}w'.format(name),common)
    print_commons('      real*8 x{0}(500),fx{0}(500,20),fx{0}tot(500)'.format(name),common)
    print_commons('      common/dist_{0}/x{0},fx{0},fx{0}tot'.format(name),common)
    print_commons('      integer m_{0},ndiv_{0} '.format(name),common)
    print_commons('      common/inp_{0}/m_{0}'.format(name),common)
    print_commons('      common/div_{0}/ndiv_{0} '.format(name),common)
    
    inc_name ='{}_dist.inc'.format(name)
    with open('inc/dists/{}'.format(inc_name),'w') as inc:
        inc.write(common.getvalue())
    inc_files.append(inc_name)
    
    print_code('C {}'.format(description[name]),binning)
    print_code('      if(m_{0}.eq.1)then'.format(name),binning)
    print_code('c {0} bin size.'.format(name),binning)
    print_code('        {0}w=({0}max-{0}min)/ndiv_{0}'.format(name),binning)
    print_code('c generate bin in {0}.'.format(name),binning)
    print_code('        do i=1,ndiv_{0}'.format(name),binning)
    print_code('          x{0}(i)={0}min+{0}w*(i-1)+{0}w/2.d0'.format(name),binning)
    print_code('        end do'.format(name),binning)
    print_code('      end if'.format(name),binning)
    
    binning_name ='{}_bin.inc'.format(name)
    with open('inc/dists/{}'.format(binning_name),'w') as binn:
        binn.write(binning.getvalue())
        
    dstring, label, symbol = words[name]
    print_commons('C {}'.format(description[name]),plotting)
    print_commons("      if(m_{0}.eq.1)then".format(name),plotting)
    print_commons("c plot distribution in {0}.".format(name),plotting)
    print_commons("        do j=1,ndiv_{0}".format(name),plotting)
    print_commons("          do i=1,it",plotting)
    print_commons("            fx{0}(j,i)=fx{0}(j,i)*avgi/cnorm(i)/{0}w".format(name),plotting)
    print_commons("          end do",plotting)
    print_commons("        end do",plotting)
    print_commons("        do j=1,ndiv_{0}".format(name),plotting)
    print_commons("          fx{0}tot(j)=0.d0 ".format(name),plotting)
    print_commons("          do i=1,it",plotting)
    print_commons("            fx{0}tot(j)=fx{0}tot(j)+fx{0}(j,i)".format(name),plotting)
    print_commons("          end do",plotting)
    print_commons("        end do",plotting)
    print_commons("        write(jh,*)'DISTRIBUTION:'",plotting)
    print_commons("        write(jh,*)'DSTRING: \"{0}\"'".format(dstring),plotting)
    print_commons("        write(jh,*)'TITLE: \"{0} Distribution\"'".format(label),plotting)
    print_commons("        write(jh,*)'AXES: \"{0} (GeV)\"',".format(symbol),plotting)
    print_commons("     &' \"d{{/Symbol s}}/{0} (pb/GeV)\"'".format(symbol),plotting)
    print_commons("        do i=1,ndiv_{0}".format(name),plotting)
    print_commons("          write(jh,964)x{0}(i),fx{0}tot(i)".format(name),plotting)
    print_commons("        end do",plotting)
    print_commons("        write(jh,*)'END'",plotting)
    print_commons("      end if",plotting)
    
    plotting_name ='{}_plot.inc'.format(name)
    with open('inc/dists/{}'.format(plotting_name),'w') as plot:
        plot.write(plotting.getvalue())
    
    print_code('      if(m_{0}.eq.1)then'.format(name),filling)
    print_code('c generate distribution in {0}.'.format(name),filling)
    print_code('        nbin=int(({0}-{0}min)/{0}w)+1'.format(name),filling)
    print_code('        if(nbin.ge.(ndiv_{0}+1))then'.format(name),filling)
    print_code('          continue',filling)
    print_code('        else if(nbin.lt.1)then',filling)
    print_code('          continue',filling)
    print_code('        else',filling)
    if name in ('eta3','eta4','cost3','cost3','y3','y4'):
        print_code('          fx{0}(nbin,it)=fx{0}(nbin,it)+hist1'.format(name),filling)
    else:
        print_code('          fx{0}(nbin,it)=fx{0}(nbin,it)+hist'.format(name),filling)
    print_code('        end if',filling)
    if name in ('eta3','eta4','cost3','cost3','y3','y4'):
        print_code('      nbin=int(({0}alt-{0}min)/{0}w)+1'.format(name),filling)
        print_code('      if(nbin.ge.(ndiv_{0}+1))then'.format(name),filling)
        print_code('        continue',filling)
        print_code('      else if(nbin.lt.1)then',filling)
        print_code('        continue',filling)
        print_code('      else',filling)
        print_code('        fx{0}(nbin,it)=fx{0}(nbin,it)+hist2'.format(name),filling)
        print_code('      end if',filling)
    print_code('      end if',filling)
    
    filling_name ='{}_fill.inc'.format(name)
    with open('inc/dists/{}'.format(filling_name),'w') as fill:
        fill.write(filling.getvalue())
    
    print_commons('      data {0}min/-999d0/,{0}max/-999d0/,ndiv_{0}/200/,m_{0}/0/'.format(name),setting)
    print 'm_{} = ndist'.format(name)
        
    
    
    
    
    
    
    
    
    
#
# with open('inc/dists/plotparams.inc','w') as sett:
#     sett.write(setting.getvalue())
with open('inc/dists.inc','w') as dists_inc:
    for inc_name in inc_files:
        dists_inc.write("      include '{}'\n".format(inc_name))
    

        
        