from coffea.nanoevents import transforms
from coffea.nanoevents.util import quote, concat
from coffea.nanoevents.schemas.base import BaseSchema, zip_forms, nest_jagged_forms


class HGCRecoSchema(BaseSchema):
  """HGCAL Reco Ntuple construction schema

  """
  def __init__(self, base_form):
    super().__init__(base_form)
    self._form["contents"] = self._build_collections(self._form["contents"])

  def _build_collections(self, branch_forms):
    collections = {
        'genpart': [],
        'pfcluster': ['hit'],
        'simcluster': ['hit'],
        'pfclusterFromMultiCl': ['hit'],
        'multiclus': ['cluster2d'],
        'vtx': [],
        'gen': ['daughter'],
        'cluster2d': ['rechits'],
        'track': ['pos'],
        'simhit': [],
        'rechit': [],
        'calopart': ['simCluster'],
        'ecalDrivenGsfele': ['pfCluster']
    }

    ## Borrowing from TreeMakerSchema, folding the ROOT::Math::Vector objects
    # into a single collection.
    composite_objects = list(
        set(k.split("/")[0] for k in branch_forms if "/" in k))
    for objname in composite_objects:
      # grab the * from "objname/objname.*"
      components = set(k[2 * len(objname) + 2:]
                       for k in branch_forms
                       if k.startswith(objname + "/"))
      if components == {
          "fCoordinates.fPt", "fCoordinates.fEta", "fCoordinates.fPhi",
          "fCoordinates.fE",
      }:
        form = zip_forms(
            {
                "pt": branch_forms.pop(f"{objname}/{objname}.fCoordinates.fPt"),
                "eta":
                branch_forms.pop(f"{objname}/{objname}.fCoordinates.fEta"),
                "phi":
                branch_forms.pop(f"{objname}/{objname}.fCoordinates.fPhi"),
                "energy":
                branch_forms.pop(f"{objname}/{objname}.fCoordinates.fE"),
            },
            objname,
            "PtEtaPhiELorentzVector",
        )
        branch_forms[objname] = form
      elif components == {
          "fCoordinates.fX", "fCoordinates.fY", "fCoordinates.fZ",
      }:
        form = zip_forms(
            {
                "x": branch_forms.pop(f"{objname}/{objname}.fCoordinates.fX"),
                "y": branch_forms.pop(f"{objname}/{objname}.fCoordinates.fY"),
                "z": branch_forms.pop(f"{objname}/{objname}.fCoordinates.fZ"),
            },
            objname,
            "ThreeVector",
        )
        branch_forms[objname] = form
      else:
        raise ValueError(f"Unrecognized class with split branches: {components}")

    def make_collection(name, items):
      if name not in branch_forms:
        branch_forms[name] = zip_forms(
            {
                k[len(name)+1:]: branch_forms.pop(k)
                for k in items
                if not k.endswith('Count')
            }, name)
        return branch_forms[name]
      else:
        collection = branch_forms[name]
        if not collection['class'].startswith('ListOffsetArray'):
          raise NotImplementedError(
              f"{name} isn't a jagged array, not sure what to do")
        for item in items:
          itemname = item[len(name) + 1:]
          collection['content']['contents'][itemname] = branch_forms.pop(
              item)['content']
        return collection

    for name in collections:
      ## Modifying collections with additional sub-structures
      subcollections = {}
      for subname in collections[name]:
        csubname = f'{name}_{subname}'
        items = sorted(
            k for k in branch_forms
            if k.startswith(csubname + '_') and not k.endswith('Count'))
        subcollections[subname] = make_collection(csubname, items)

      items = sorted(k for k in branch_forms if k.startswith(name + "_"))
      # string the subcollections stuff
      for subname in collections[name]:
        items = [
            k for k in items
            if not k.startswith(f"{name}_{subname}") or k.endswith('Count')
        ]

      make_collection(name, items)

      for subname in collections[name]:
        nest_jagged_forms(
            branch_forms[name],
            branch_forms.pop(f"{name}_{subname}"),
            f"{subname}Count",
            f"{subname}",
        )


# Adding more stuff to the return output

    return branch_forms

  @property
  def behavior(self):
    """Behaviors necessary to implement this schema"""
    from coffea.nanoevents.methods import base, vector

    behavior = {}
    behavior.update(base.behavior)
    behavior.update(vector.behavior)
    return behavior
