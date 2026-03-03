from extractors.decimer_extractor import DECIMERExtractor


class _StubPredictor:
    def __init__(self, output):
        self.output = output
        self.calls = []

    def __call__(self, image_path):
        self.calls.append(image_path)
        return self.output


def _make_extractor(output):
    extractor = DECIMERExtractor()
    extractor.predictor = _StubPredictor(output)
    extractor._loaded = True
    return extractor


def test_decimer_extract_preserves_dict_confidence():
    predictor_output = {
        "smiles": "CN1CCC2=CC=CC=C2C1",
        "confidence": 0.83,
    }

    extractor = _make_extractor(predictor_output)
    result = extractor.extract("/tmp/fake-decimer-image.png")

    assert result["success"] is True
    assert result["smiles"] == predictor_output["smiles"]
    assert result["confidence"] == predictor_output["confidence"]


def test_decimer_extract_preserves_tuple_confidence():
    smiles = "CN1CCC2=CC=CC=C2C1"
    extractor = _make_extractor((smiles, 0.91))
    result = extractor.extract("/tmp/fake-decimer-image.png")

    assert result["success"] is True
    assert result["smiles"] == smiles
    assert result["confidence"] == 0.91


def test_decimer_extract_falls_back_when_confidence_missing():
    smiles = "C"

    extractor = _make_extractor(smiles)
    result = extractor.extract("/tmp/fake-decimer-image.png")

    assert result["success"] is True
    assert result["smiles"] == smiles
    assert result["confidence"] == extractor._estimate_confidence(smiles)
