import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import AsyncMock, patch

from phospho_group_scraper import load_scrape_state, resolve_protein_name_on_page, run_site_batch


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


class ScrapeCheckpointTests(unittest.IsolatedAsyncioTestCase):
    async def test_completed_site_ids_are_skipped_when_resuming_resolved_id_scrape(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_path = Path(tmpdir) / "scrape_state.json"
            state_path.write_text(
                json.dumps(
                    {
                        "proteins": {
                            "570": {
                                "protein_id": 570,
                                "status": "partial",
                                "discovered_site_ids": [111, 222],
                                "completed_site_ids": [111],
                                "failed_site_ids": [],
                                "cloudflare_blocked_site_ids": [],
                                "site_outputs": {"111": "existing.csv"},
                                "site_errors": {},
                            }
                        }
                    }
                ),
                encoding="utf-8",
            )
            state = load_scrape_state(state_path)

            with patch("phospho_group_scraper.main", new=AsyncMock(return_value="new.csv")) as fake_main:
                result = await run_site_batch(
                    [111, 222],
                    delay=0,
                    continue_on_error=True,
                    protein_id=570,
                    scrape_state=state,
                    scrape_state_path=state_path,
                )

            fake_main.assert_awaited_once_with(222)
            self.assertEqual(result["status"], "complete")
            saved = load_scrape_state(state_path)
            protein_state = saved["proteins"]["570"]
            self.assertEqual(protein_state["completed_site_ids"], [111, 222])
            self.assertEqual(protein_state["site_outputs"]["222"], "new.csv")


if __name__ == "__main__":
    unittest.main()
