import unittest

from phospho_group_scraper import resolve_protein_name_on_page


class FakeLookupPage:
    def __init__(self):
        self.url = ""

    async def goto(self, url, **_kwargs):
        self.url = url
        if "proteinAction.action?id=570" in url:
            self.url = "https://www.phosphosite.org/confirmAction"

    async def title(self):
        return "PhosphoSitePlus"

    def locator(self, _selector):
        return self

    async def inner_text(self, timeout=None):
        return "PhosphoSitePlus"

    async def wait_for_timeout(self, _timeout):
        return None

    async def wait_for_selector(self, _selector, timeout=None):
        return None

    async def eval_on_selector_all(self, selector, _script):
        if selector != "#simpleSearchResultsTable tbody tr":
            return []
        return [
            {
                "protein": "Akt1",
                "gene": "AKT1",
                "organism": "human",
                "links": [
                    {
                        "href": "https://www.phosphosite.org/proteinAction.action?id=570&showAllSites=true",
                        "text": "Akt1",
                    }
                ],
            }
        ]


class ProteinLookupTests(unittest.IsolatedAsyncioTestCase):
    async def test_human_candidate_id_is_saved_when_validation_redirects_to_confirm_action(self):
        result = await resolve_protein_name_on_page(FakeLookupPage(), "AKT1", "human")

        self.assertEqual(result["protein_id"], 570)
        self.assertEqual(result["organism"], "human")
        self.assertEqual(result["url"], "https://www.phosphosite.org/proteinAction.action?id=570&showAllSites=true")
        self.assertEqual(result["source"], "human_search_result_url_confirm_redirect")


if __name__ == "__main__":
    unittest.main()
