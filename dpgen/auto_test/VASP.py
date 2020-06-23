import os
import json
from dpgen import dlog
from dpgen.util import sepline
import dpgen.auto_test.lib.vasp as vasp
from dpgen.auto_test.Task import Task
from dpgen.generator.lib.vasp import incar_upper
from dpdata import LabeledSystem
from pymatgen.io.vasp import Incar, Kpoints


class VASP(Task):
    def __init__(self,
                 inter_parameter,
                 path_to_poscar):
        self.inter = inter_parameter
        self.incar = inter_parameter['incar']
        self.potcars = inter_parameter['potcars']
        default_potcar_prefix = ''
        self.potcar_prefix = inter_parameter.get('potcar_prefix',default_potcar_prefix)
        self.path_to_poscar = path_to_poscar

    def make_potential_files(self,
                             output_dir):
        with open(os.path.join(output_dir, 'POTCAR'), 'w') as fp:
            for ii in self.potcars:
                with open(os.path.join(self.potcar_prefix, self.potcars[ii]), 'r') as fin:
                    for line in fin:
                        print(line.strip('\n'), file=fp)

        with open(os.path.join(output_dir, 'inter.json'), 'w') as fp:
            json.dump(self.inter, fp, indent=4)

    def make_input_file(self,
                        output_dir,
                        task_type,
                        task_param):
        sepline(ch=output_dir)
        with open(os.path.join(output_dir, 'task.json'), 'w') as fp:
            json.dump(task_param, fp, indent=4)

        assert (os.path.exists(self.incar)), 'no INCAR file for relaxation'
        relax_incar_path = os.path.abspath(self.incar)
        incar = incar_upper(Incar.from_file(relax_incar_path))

        if 'ISIF' in incar:
            isif = incar.get('ISIF')
        else:
            isif = 3

        if 'NSW' in incar:
            nsw = incar.get('NSW')
        else:
            nsw = 200

        try:
            kspacing = incar.get('KSPACING')
        except KeyError:
            raise RuntimeError("KSPACING must be given in INCAR")

        if 'KGAMMA' in incar:
            kgamma = incar.get('KGAMMA')
        else:
            kgamma = False


        if task_type in ['relaxation', 'vacancy', 'interstitial']:
            isif = 3

        if task_type == 'eos':
            if 'change_box' in task_param and not task_param['change_box']:
                isif = 2
            else:
                isif = 4

        if task_type == 'elastic':
            isif = 2

        if task_type == 'surface':
            if 'static-opt' in task_param and task_param['static-opt']:
                nsw = 0
            elif 'change_box' in task_param and task_param['change_box']:
                isif = 4
            else:
                isif = 2

        if task_type == 'static' \
                or ('reprod_opt' in task_param and task_param['reprod_opt']):
            nsw = 0

        if not ('ISIF' in incar and incar.get('ISIF') == isif):
            dlog.info("%s setting ISIF to %d" % ( self.make_input_file.__name__, isif))
            incar['ISIF'] = isif

        if not ('NSW' in incar and incar.get('NSW') == nsw):
            dlog.info("%s setting NSW to %d" % ( self.make_input_file.__name__, nsw))
            incar['NSW'] = nsw


        if 'ediff' in task_param:
            dlog.info("%s setting ediff to %s" % (self.make_input_file.__name__, task_param['ediff']))
            incar['EDIFF'] = task_param['ediff']

        if 'ediffg' in task_param:
            dlog.info("%s setting ediffg to %s" % ( self.make_input_file.__name__, task_param['ediffg']))
            incar['EDIFFG'] = task_param['ediffg']

        if 'encut' in task_param:
            dlog.info("%s setting encut to %s" % (self.make_input_file.__name__, task_param['encut']))
            incar['ENCUT'] = task_param['encut']


        if 'kspacing' in task_param:
            dlog.info("%s setting kspacing to %s" %( self.make_input_file.__name__, task_param['kspacing']))
            incar['KSPACING'] = task_param['kspacing']

        if 'kgamma' in task_param:
            dlog.info("%s setting kgamma to %s" % (  self.make_input_file.__name__, task_param['kgamma']))
            incar['KGAMMA'] = task_param['kgamma']

        incar.write_file(os.path.join(output_dir,'INCAR'))
        ret = vasp.make_kspacing_kpoints(self.path_to_poscar, kspacing, kgamma)
        kp = Kpoints.from_string(ret)
        kp.write_file(os.path.join(output_dir, "KPOINTS"))

    def compute(self,
                output_dir,inter_param=None):
        outcar = os.path.join(output_dir, 'OUTCAR')
        if not os.path.isfile(outcar):
            dlog.warning("cannot find OUTCAR in " + output_dir + " skip")
            return None
        else:
            ls=LabeledSystem(outcar)
            if len(ls)>0:
               force  = ls.sub_system([-1]).data['forces'][0].tolist()
               energy = ls.sub_system([-1]).data['energies'][0].tolist()
               virials= ls.sub_system([-1]).data['virials'][0].tolist()
               return {"energy": energy, "force": force, "virials":virials}

    def forward_files(self):
        return ['INCAR', 'POSCAR', 'POTCAR']

    def forward_common_files(self):
        return ['INCAR', 'POTCAR']

    def backward_files(self):
        return ['OUTCAR', 'outlog', 'CONTCAR', 'OSZICAR']