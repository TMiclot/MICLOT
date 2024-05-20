#!/usr/bin/env python3

"""
  ___      ___   __     ______   ___        ______  ___________  
 |"  \    /"  | |" \   /" _  "\ |"  |      /    " \("     _   ") 
  \   \  //   | ||  | (: ( \___)||  |     // ____  \)__/  \\__/  
  /\   \/.    | |:  |  \/ \     |:  |    /  /    ) :)  \\_ /     
 |: \.        | |.  |  //  \ _   \  |___(: (____/ //   |.  |     
 |.  \    /:  | /\  |\(:   _) \ ( \_|:  \\        /    \:  |     
 |___|\__/|___|(__\_|_)\_______) \_______)\"_____/      \__|

 Molecular InteraCtion anaLysis tOolkiT
 _______________________________________
 
 A collection of tools to analyse protein-protein interactions.
"""

__author__ = 'Tom MICLOT  <tom.miclot@jh-inst.cas.cz>'
__license__ = "xxxx"
__version__ = "Version: 1.0 -- jj/mm/2024"

__all__ = ['database', 'complex_binding', 'utilities', 'cys_bridges', 'interactions', 'coulomb_lj']


#=====| Import all scripts|=====
from . import database
from . import utilities
from . import complex_binding
from . import cys_bridges
from . import interactions
from . import coulomb_lj



#=====| Module end |=====
if __name__ == "__main__":
    main()

