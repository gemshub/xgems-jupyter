"""
A simple object oriented interface to control gems in python with pythonic 
naming convention and dictionaries as a more convenient wrapper on pyxGEMS.
Created by R.A.Patel (c) 2021.

This will be superseesed in 2025 by new ChemicalEngingeDicts a dictionaries based ChemicalEngine
"""
#%%import libraries
import xgems
import numpy as np
#%%
class GEMS(object):
    """
    Gems interface in calculator format for easy use
    """
    def __init__(self,input_file, reset_calc=True, cold_start=True):
        """
        Initialization of gems calculator
        reset_calc: True will reset b vector to zero
        """
        self.input_file = input_file
        self.gem = xgems.ChemicalEngine(input_file)
        self.T = self.gem.temperature()
        self.P = self.gem.pressure()
        self.b = np.array(self.gem.elementAmounts())
        if cold_start: self.cold_start()
        self.equilibrate()
        self.nelements = self.gem.numElements()
        self.nphases = self.gem.numPhases()
        self.nspecies = self.gem.numSpecies()
        self.element_names = []
        self.element_molar_masses ={}
        elemolarmass = self.gem.elementMolarMasses()
        for i in range(self.nelements):
            self.element_names.append(self.gem.elementName(i))
            self.element_molar_masses[self.gem.elementName(i)] = elemolarmass[i]
        self.species_names =[]
        for i in range(self.nspecies):
            self.species_names.append(self.gem.speciesName(i))
        self.phase_names = []
        for i in range(self.nphases):
            self.phase_names.append(self.gem.phaseName(i))
        #dictionary containing species in phase
        self.species_in_phase = {}
        for i in range(self.nphases):
            if self.gem.numSpeciesInPhase(i) > 0: 
                  self.species_in_phase[self.phase_names[i]] = self.species_names[self.gem.indexFirstSpeciesInPhase(i):
                      self.gem.indexFirstSpeciesInPhase(i)+self.gem.numSpeciesInPhase(i)]
        molar_mass = self.gem.speciesMolarMasses()
        self.species_molar_mass={}
        for i in range(self.nspecies):
            self.species_molar_mass[self.species_names[i]] = molar_mass[i]
        self.species_molar_volumes = {}
        for i in range(self.nspecies):
            self.species_molar_volumes[self.species_names[i]] = self.gem.standardMolarVolume(i) 
        self.formulaMatrix=self.gem.formulaMatrix().T
        if reset_calc:self.clear()

    def clear(self):
        self.b[:]=1e-12
    
    @property
    def aq_composition(self):
        """
        aq solution composition in mole/kgW
        """
        out = {}
        molalities =  self.gem.speciesMolalities()
        aq_species_names=self.species_in_phase['aq_gen'][:-1]
        for name in aq_species_names:
            out[name] =molalities[self.species_names.index(name)]
        return out
    
    def reset_aq_composition(self):
        """
        Removes all ions from aqueous solution
        be careful as this will also remove water i.e H+ and OH-
        """
        self.b -= self.gem.elementAmountsInPhase(self.phase_names.index('aq_gen'))
    
    @property
    def phase_amounts(self):
        """
        return phase amounts in moles
        """
        pamt = self.gem.phaseAmounts()
        out = {}      
        for name in self.phase_names:
            out[name]=pamt[self.phase_names.index(name)]
        return out
    
    @property
    def species_amounts(self):
        """
        returns species amounts in moles
        """
        spamt = self.gem.speciesAmounts()
        out = {}      
        for name in self.species_names:
            out[name]=spamt[self.species_names.index(name)]
        return out


    @property
    def solid_mass_frac(self):
        """
        m/m ratio for  solid phases 
        """
        mfrac = self.gem.phaseMasses()/np.sum(self.gem.phaseMasses())
        out = {}      
        for name in self.phase_names:
            out[name]=mfrac[self.phase_names.index(name)]
        return out
    
    @property
    def solid_volume_frac(self):
        """
        m3/m3 ratio for solid phases
        """
        out = self.phase_volume_frac
        del out["aq_gen"],out["gas_gen"]
        return out
    
    @property
    def aq_volume_frac(self):
        """
        Volume fraction of water in the system
        """
        return self.phase_volume_frac["aq_gen"]
    
    @property
    def phase_volumes(self):
        """
        returns absolute volume of phases in the system
        """
        v = self.gem.phaseVolumes()
        out = {}      
        for name in self.phase_names:
            out[name]=v[self.phase_names.index(name)]
        return out
    
    @property
    def phase_masses(self):
        mass = self.gem.phaseMasses()
        out = {}      
        for name in self.phase_names:
            out[name]=mass[self.phase_names.index(name)]
        return out
    
    @property
    def phase_volume_frac(self):
        """
        volume fraction of all phases in the system
        """
        vfrac = self.gem.phaseVolumes()/self.system_volume
        out = {}      
        for name in self.phase_names:
            out[name]=vfrac[self.phase_names.index(name)]
        return out
    
    def equilibrate(self):
        """
        runs equilibriation of the current recipe
        """
        outcode= self.gem.equilibrate(self.T,self.P,self.b)
        return self._status_encoder[outcode]

    def cold_start(self):
        self.gem.setColdStart()
    
    def warm_start(self):
        self.gem.setWarmStart()
        
    def add_multiple_species_amt(self,input_dict,units="moles"):
        """
        add species amount in the system useful for adding aqueous solution composition
        units= moles,kg,m3
        """
        for name, val in input_dict.items():
            self.add_species_amt(name,val,units)

    def add_species_amt(self,species,val,units="moles"):
        """
        add species amount in the system useful for adding aqueous solution composition
        units= moles,kg,m3
        """
        species_idx = self.species_names.index(species)
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        self.b += self.formulaMatrix[species_idx]*val
    
    def add_element_amt(self,element_name,val, units = "moles"):
        """
        add element amount in the system
        units = moles, kg
        """
        if units  == "kg":
            val /= self.element_molar_masses[element_name]  
        self.b[self.element_names.index(element_name)]+=val
        
    def add_amt_from_formula(self,formula, val, units="moles"):
        """
        add element amount using user defined formula
        units = moles,kg
        """
        if units == "kg":
            molarmass =0 
            for element in formula.keys():
                molarmass += formula[element] * self.element_molar_masses[element]
            val/=molarmass
        for element in formula.keys():
            self.add_element_amt(element, val * formula[element])
    
    def multiple_species_lower_bound(self, input_dict,units="moles"):
        """
        constrain species amount to a specified lower bound 
        units= moles,kg,m3
        """
        for species,val in input_dict.items():
            self.species_lower_bound(species,val,units)
        
    def multiple_species_upper_bound(self, input_dict,units="moles"):
        """
        constrain species amount to a specified lower bound 
        units= moles,kg,m3
        """
        for species,val in input_dict.items():
            self.species_upper_bound(species,val,units)
    
    def species_lower_bound(self, species,val,units="moles"):
        """
        constrain species amount to a specified lower bound 
        units= moles,kg,m3
        """
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        species_idx = self.species_names.index(species)
        self.gem.setSpeciesLowerLimit(species_idx,val)
    
    def species_upper_bound(self,species,val,units="moles"):
        """
        constrain species amount to a specified upper bound
        units= moles,kg,m3
        """
        if units == "kg":
            val/=self.species_molar_mass[species]
        if units == "m3":
            val/=self.species_molar_volumes[species]
        self.gem.setSpeciesUpperLimit(species,val)
        
    def supress_phase(self,phase_name):
        """
        supresses phase in calculation
        """
        for species in self.species_in_phase[phase_name]:
            self.supress_species(species)
    
    def supress_multiple_phases(self,phase_name_list):
        """
        supresses multiple phase in calculation given in list
        """
        for phase in phase_name_list:
            self.supress_phase(phase)
    
    def supress_multiple_species(self,species_list):
        """
        supresses multiple species in calculation given in list
        """
        for species in species_list:
            self.supress_species(species)
    
    def supress_species(self,species_name):
        """
        supresses species in calculation
        """
        self.species_lower_bound(species_name,0)
        self.species_upper_bound(species_name,1e-15)
        
    def activate_phase(self,phase_name):
        """
        activate supressed phase
        """
        for species in self.species_in_phase[phase_name]:
            self.activate_species(species)
    
    def activate_multiple_phases(self,phase_name_list):
        """
        activate multiple supressed phases given in list
        """
        for phase in phase_name_list:
            self.activate_phase(phase)
    
    def activate_multiple_species(self,species_list):
        """
        activate multiple supressed species given in the list
        """
        for species in species_list:
            self.supress_species(species)
    
    def activate_species(self,species_name):
        """
        activate a supressed species in phase
        """
        self.species_lower_bound(species_name,0)
        self.species_upper_bound(species_name,1e6)
            
    @property
    def pH(self):
        """
        returns pH of the solution
        """
        return self.gem.pH()
        
    @property
    def pE(self):
        """
        returns pE of the solution
        """
        return self.gem.pe()
        
    @property
    def ionic_strength(self):
        """
        returns ionic strength of the solution
        """
        return self.gem.ionicStrength()
    
    @property
    def system_volume(self):
        """
        returns volume of the system
        """
        return self.gem.systemVolume()

    @property
    def system_mass(self):
        """
        returns mass of the system
        """
        return self.gem.systemMass()
    
    @property
    def phase_molar_volume(self):
        """
        returns molar volume of phases
        """
        phase_mvol = {}
        for i in range(self.nphases):
            phase_mvol[self.phase_names[i]]=self.gem.phaseMolarVolume(i)
        return phase_mvol
    
    @property
    def phase_sat_indices(self):
        """
        returns saturation indices of phases
        """
        phaseSatIndices= self.gem.phaseSatIndices()
        satIndices = {}
        for i in range(self.nphases):
            satIndices[self.phase_names[i]]=phaseSatIndices[i]
        return satIndices
    
    @property
    def _status_encoder(self):
        out = {}
        out[0]="No GEM re-calculation needed"
        out[1]= "Need GEM calculation with LPP (automatic) initial approximation (AIA)"
        out[2]="OK after GEM calculation with LPP AIA"
        out[3]="Bad (not fully trustful) result after GEM calculation with LPP AIA"
        out[4]="Failure (no result) in GEM calculation with LPP AIA"
        out[5]="Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation (full DATABR lists only)"
        out[6]="OK after GEM calculation with SIA"
        out[7]="Bad (not fully trustful) result after GEM calculation with SIA"
        out[8]="Failure (no result) in GEM calculation with SIA"
        out[9]="Terminal error has occurred in GEMS3K (e.g. memory corruption). Restart is required."
        return out         
    
if __name__ == "__main__":    
    """
    Dummy cement hydration case
    """
    def plot(results):
        import matplotlib.pylab as plt
        prev = 0
        ncat=len((np.nonzero(np.fromiter(results.values(),dtype=float)))[0])
        category_colors = plt.get_cmap('tab10')(np.linspace(0., 1, ncat+1))
        i = 0
        for name,val in results.items():
            if val !=0:
                i+=1
                plt.barh(0, val , left=prev, height=0.2,
                        label=name,color =category_colors[i])
                prev += val
        plt.legend(ncol=3,bbox_to_anchor=(0, 2.2), loc='upper left')
        ax = plt.gca()
        ax.axes.yaxis.set_visible(False)
        plt.axis("image")
        plt.show()
 
    DoH =0.2  #assumed for all phases for test case here
    C3S=64.6  #g/100g cement
    C2S=9.3   #g/100g cement
    C3A=7.4   #g/100g cement
    C4AF=7.8  #g/100g cement  
    wc= 0.5   #g/100g cement
    input_file = 'gems_files/CemHyds-dat.lst'
    gem = GEMS(input_file)
    species_dict = {"C3S":C3S *1e-3 ,
                    "C2S":C2S*1e-3,
                    "C3A":C3A*1e-3 ,
                    "C4AF":C4AF*1e-3,
                    }
    gem.add_amt_from_formula({"H":2,"O":1},wc * 100 * 1e-3,units="kg")
    gem.add_species_amt("O2",1e-6)
    gem.add_multiple_species_amt(species_dict,units = "kg")
#    set lower bound to limit dissolution for unreacted phases
    for name,val in species_dict.items():
        if name!="H2O@":
            gem.species_lower_bound(name,val*(1-DoH),units="kg")
    print(gem.equilibrate())
    plot(gem.phase_volume_frac)

   

     