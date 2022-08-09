
import torch
from bgmol.datasets import WaterDimerFlexibleTIP3P


def test_torch_datasets(tmpdir):
    dataset = WaterDimerFlexibleTIP3P(root=str(tmpdir))
    trainset, valset, testset = dataset.torch_datasets(val_fraction=0.4, fields=["xyz", "energies"])
    assert len(trainset) == 30_000
    assert len(valset) == 20_000
    assert len(testset) == 0
    xyz, energies = trainset[0]
    assert xyz.shape == torch.Size([6, 3])
    assert energies.shape == torch.Size([])
    trainset2, valset2, testset2 = dataset.torch_datasets(random_seed=1, val_fraction=0.4, fields=["xyz", "energies"])
    trainset3, valset3, testset3 = dataset.torch_datasets(random_seed=2, val_fraction=0.4, fields=["xyz", "energies"])
    assert torch.allclose(trainset[0][0], trainset2[0][0])
    assert torch.allclose(trainset[0][1], trainset2[0][1])
    assert not torch.allclose(trainset[0][0], trainset3[0][0])
